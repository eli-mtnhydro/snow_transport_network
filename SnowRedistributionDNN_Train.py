# Copyright Eli Boardman 2024
# All rights reserved

import os
import torch
import pandas as pd
import math
import csv

os.chdir('C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/')

device = torch.device('cuda')

#%%

class RedistributionNetwork(torch.nn.Module):
    def __init__(self, nX, nY, initial_storage_data, input_scale):
        super(RedistributionNetwork, self).__init__()
        
        self.efficiency = torch.nn.Parameter(torch.tensor(1.0).to(device))
        self.input_scale = torch.tensor(input_scale).to(device)
        
        self.nX = nX
        self.nY = nY
        
        # Pad nY with 0s to make weights matrix a multiple of 8 for tensor cores
        if (nY % 8) > 0:
            self.nYpad = self.nY - (self.nY % 8) + 8
        else:
            self.nYpad = self.nY
        
        self.initial_flux =  torch.zeros(self.nYpad).to(device)
        self.initial_storage = torch.zeros(self.nYpad, nX).to(device)
        
        self.direction_distances = []
        self.connection_mask = []
        
        self.layers = torch.nn.ModuleList()
        
        for x in range(nX):
            layer = torch.nn.Linear(self.nYpad, self.nYpad, bias=False)
            layer.weight.data.fill_(0)
            
            direction_distance = torch.zeros_like(layer.weight.data).to(device)
            connection_mask = torch.zeros_like(layer.weight.data).to(device)
            
            for y in range(nY):
                self.initial_storage[y, x] = initial_storage_data[y, x]
                
                # All weights remain zero in last layer to enforce no-transport boundary condition
                if x < (nX - 1):
                    if y > 0:
                        layer.weight.data[y-1, y] = 0.01
                        connection_mask[y-1, y] = 1
                        direction_distance[y-1, y] = math.sqrt(2)
                    layer.weight.data[y, y] = 0.1
                    connection_mask[y, y] = 1
                    direction_distance[y, y] = 1
                    if y < nY-1:
                        layer.weight.data[y+1, y] = 0.01
                        connection_mask[y+1, y] = 1
                        direction_distance[y+1, y] = math.sqrt(2)
            
            self.layers.append(layer)
            self.direction_distances.append(direction_distance)
            self.connection_mask.append(connection_mask)
        
        self.input_sum = self.initial_storage.sum()
        
    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
        return x
    
    # For training
    def compute_storage_and_cost(self):
        total_loss = torch.tensor(0.0,requires_grad=True).to(device)
        final_storage = torch.zeros_like(self.initial_storage).to(device)
        flux_cost = torch.tensor(0.0,requires_grad=True).to(device)
        incoming_flux = self.initial_flux
        
        for x, layer in enumerate(self.layers):
            
            # Available mass from incoming flux and local initial storage
            total_avail = incoming_flux + (self.initial_storage[:,x] * self.input_scale)
            
            # Remaining (final) mass storage
            final_storage[:,x] = total_avail * (1 - layer.weight.sum(dim=1))
            
            # Determine incoming flux to next layer
            incoming_flux = total_avail @ layer.weight
            
            # Account for loss
            total_loss = total_loss + sum(incoming_flux * (1.0 - self.efficiency))
            incoming_flux = incoming_flux * self.efficiency
            
            flux_cost = flux_cost + sum(total_avail @ (layer.weight * self.direction_distances[x]))
        
        loss_frac = total_loss / (self.input_sum * self.input_scale)
        return final_storage[0:self.nY, :], flux_cost, loss_frac
    
    # For creating final results
    def compute_storage_and_flux(self):
        net_loss = torch.zeros_like(self.initial_storage).to(device)
        final_storage = torch.zeros_like(self.initial_storage).to(device)
        total_flux = torch.zeros_like(self.initial_storage).to(device)
        incoming_flux = self.initial_flux
        
        for x, layer in enumerate(self.layers):
            
            # Available mass from incoming flux and local initial storage
            total_avail = incoming_flux + (self.initial_storage[:,x] * self.input_scale)
            
            # Remaining (final) mass storage
            final_storage[:,x] = total_avail * (1 - layer.weight.sum(dim=1))
            
            # Determine incoming flux to next layer
            incoming_flux = total_avail @ layer.weight
            
            # Account for loss
            net_loss[:,x] = incoming_flux * (1.0 - self.efficiency)
            incoming_flux = incoming_flux * self.efficiency
            
            total_flux[:,x] = incoming_flux
            
        return final_storage[0:self.nY, :], total_flux[0:self.nY, :], net_loss[0:self.nY, :]
    
#%%

class WeightClipper(object):
    def __init__(self, max_angle_scale, nY, nYpad, min_weight, max_weight, leaky_factor):
        self.max_angle_scale = max_angle_scale
        self.nY = nY
        
        self.upper_diagonal_mask = torch.zeros(nYpad, nYpad, dtype=torch.bool).to(device)
        self.lower_diagonal_mask = torch.zeros(nYpad, nYpad, dtype=torch.bool).to(device)
        
        self.upper_diagonal_mask[(torch.arange(0, nY - 1), torch.arange(1, nY))] = True
        self.lower_diagonal_mask[(torch.arange(1, nY), torch.arange(0, nY - 1))] = True
        
        self.min_weight = min_weight
        self.max_weight = max_weight
        self.leaky_factor = leaky_factor
        
    def __call__(self, module):
        if hasattr(module, 'weight'):
            weights = module.weight.data
            
            # Clip to max allowed range
            weights = weights.clamp(self.min_weight,self.max_weight)
            
            # Normalize the weights of off-diagonal neurons to respect maximum deflection angle
            diagonal_weights = weights.diagonal()
            upper_diagonal_weights = (self.upper_diagonal_mask * weights).sum(dim=1)
            lower_diagonal_weights = (self.lower_diagonal_mask * weights).sum(dim=1)
            upper_diagonal_weights[upper_diagonal_weights < lower_diagonal_weights] = self.leaky_factor
            lower_diagonal_weights[upper_diagonal_weights > lower_diagonal_weights] = self.leaky_factor
            upper_diagonal_weights = upper_diagonal_weights.clamp(max=self.max_angle_scale * diagonal_weights)
            lower_diagonal_weights = lower_diagonal_weights.clamp(max=self.max_angle_scale * diagonal_weights)
            weights[self.upper_diagonal_mask] = upper_diagonal_weights[0:(self.nY-1)]
            weights[self.lower_diagonal_mask] = lower_diagonal_weights[1:self.nY]
            
            # Normalize the weights for neurons that exceed a sum of 100%
            exceed_mask = weights.sum(dim=1) > self.max_weight
            if exceed_mask.any():
                active_weights = weights[exceed_mask]
                normalized_weights = self.max_weight * active_weights / active_weights.sum(dim=1, keepdim=True)
                weights[exceed_mask] = normalized_weights
            
            module.weight.data = weights
        
        if hasattr(module, 'efficiency'):
            module.efficiency.data = module.efficiency.data.clamp(min=0.0,max=1.0)

#%%

def load_data(file_path):
    data = pd.read_csv(file_path, header=None)
    return torch.tensor(data.values, dtype=torch.float32)

#%%

def main():
    
    # User inputs
    scenario_name = 'Smoothed2x21_22.5deg_0.1loss'
    n_epochs = 100000
    max_deflection_angle = 22.5 # Degrees up or down relative to across
    target_loss_frac = 0.1 # Fraction of total accumulation that is sublimated
    target_loss_tol = 0.001 # Must be within this percentage of target_loss_frac
    target_rmse = 0.08
    max_cell_err_multiple = 2
    max_cell_err_multiple_buffer = max_cell_err_multiple * 0.99 # To avoid continually re-training cells close to threshold
    nSkip_for_max_err = 10
    target_bias = 0.01
    depth_thresh_spacing = 0.1
    min_weight = 0
    max_weight = 0.99
    leaky_factor = 1e-4 # To keep nonzero gradients on up/down tridiagonal weights
    
    input_scale = 1 / (1 - target_loss_frac)
    
    # Compute max deflection between diagonal and off-diagonal weights
    tan_max_deflection_angle = math.tan(max_deflection_angle * math.pi/180)
    max_angle_scale = tan_max_deflection_angle / (1 - tan_max_deflection_angle)
    
    # Load data and create target tensor
    target_data = load_data('MassDistribution_Reference.csv').to(device)
    initial_storage_data = load_data(f'{scenario_name}/MassDistribution_Counterfactual.csv').to(device)
    
    nX = len(target_data[0]) # Number of layers
    nY = len(target_data)    # Number of neurons per layer
    nTotal = nX * nY
    
    max_err_layer_mask = torch.ones(nY, nX).to(device)
    max_err_layer_mask[:,(nX-1)] = 0 # Last layer is not able to get rid of excess initial mass
    if nSkip_for_max_err > 0:
        for x in range(nSkip_for_max_err):
            max_err_layer_mask[:,x] = 0 # Early layers may not have enough mass available for import
    
    model = RedistributionNetwork(nX, nY, initial_storage_data, input_scale).to(device)
    clipper = WeightClipper(max_angle_scale, nY, model.nYpad, min_weight, max_weight, leaky_factor)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=0.0)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5,
                                                           patience=1000, threshold=1e-3, threshold_mode="rel",
                                                           cooldown=0, min_lr=1e-5)
    scaler = torch.cuda.amp.GradScaler()
    
    # Track training performance
    lr_list = []
    rmse_list = []
    rmse_binned_list = []
    bias_list = []
    flux_cost_list = []
    
    # Training loop
    torch.manual_seed(42)
    for epoch in range(n_epochs):
        
        with torch.autocast(device_type='cuda', dtype=torch.float16):
            bias_val = torch.tensor(0.0,requires_grad=True).to(device)
            binned_rmse_combo = torch.tensor(0.0,requires_grad=True).to(device)
            final_storage, flux_cost, loss_frac = model.compute_storage_and_cost()
            flux_cost = flux_cost / nTotal
            storage_errs = target_data - final_storage
            
            rmse_val = torch.sqrt(torch.mean(torch.square(storage_errs)))
            
            err_mask = torch.abs(storage_errs * max_err_layer_mask) > (target_rmse * max_cell_err_multiple)
            if err_mask.any():
                err_mask_extra = torch.abs(storage_errs * max_err_layer_mask) > (target_rmse * max_cell_err_multiple_buffer)
                binned_rmse_combo = torch.sqrt(torch.mean(torch.square(storage_errs[err_mask_extra])))
            
            max_depth = torch.max(final_storage).item()
            depth_thresh = 0
            nBiasDepths = 0
            while depth_thresh < max_depth:
                depth_thresh_mask = (final_storage > depth_thresh) * (final_storage < (depth_thresh + depth_thresh_spacing))
                depth_thresh += depth_thresh_spacing
                if depth_thresh_mask.any():
                    binned_depth_errs = target_data[depth_thresh_mask] - final_storage[depth_thresh_mask]
                    bias_val = bias_val + torch.abs(torch.mean(binned_depth_errs))
                    nBiasDepths += 1
            bias_val = bias_val / nBiasDepths
        
        target_loss_err = torch.abs(target_loss_frac - loss_frac) / target_loss_frac
        
        # Multi-phase loss function
        if rmse_val < target_rmse and binned_rmse_combo < target_rmse and bias_val < target_bias and target_loss_err < target_loss_tol:
            epoch_loss = flux_cost
        else:
            if rmse_val < target_rmse:
                if binned_rmse_combo > 0:
                    epoch_loss = binned_rmse_combo
                    if bias_val > target_bias:
                        epoch_loss = epoch_loss + bias_val * 0.01 * (epoch_loss.item() / bias_val.item())
                    if target_loss_err > target_loss_tol:
                        epoch_loss = epoch_loss + target_loss_err * 0.01 * (epoch_loss.item() / target_loss_err.item())
                elif bias_val > target_bias:
                    epoch_loss = bias_val
                    if target_loss_err > target_loss_tol:
                        epoch_loss = epoch_loss + target_loss_err * 0.01 * (epoch_loss.item() / target_loss_err.item())
                else:
                    epoch_loss = target_loss_err
            else:
                epoch_loss = rmse_val
                if binned_rmse_combo > 0:
                    epoch_loss = epoch_loss + binned_rmse_combo * 0.01 * (epoch_loss.item() / binned_rmse_combo.item())
                if bias_val > target_bias:
                    epoch_loss = epoch_loss + bias_val * 0.01 * (epoch_loss.item() / bias_val.item())
                if target_loss_err > target_loss_tol:
                    epoch_loss = epoch_loss + target_loss_err * 0.01 * (epoch_loss.item() / target_loss_err.item())
            # Add small amount of flux_cost gradient
            if flux_cost > 0:
                epoch_loss = epoch_loss + flux_cost * 0.01 * (epoch_loss.item() / flux_cost.item())
        
        scaler.scale(epoch_loss).backward()
        
        # Zero out gradients for non-connected weights (initial weight == 0)
        for x in range(nX):
            model.layers[x].weight.grad *= model.connection_mask[x]
        
        scaler.step(optimizer)
        scale = scaler.get_scale()
        scaler.update()
        if not (scale > scaler.get_scale()):
            scheduler.step(epoch_loss)
        optimizer.zero_grad()
        model.apply(clipper)
        
        current_lr = optimizer.param_groups[0]['lr']
        if epoch < 10 or (epoch < 100 and epoch % 10 == 0) or (epoch < 1000 and epoch % 100 == 0) or epoch % 1000 == 0:
            print(f'Epoch {epoch}, Lr: {current_lr:.3g}, RMSE: {rmse_val.item():.3g}, Bias: {bias_val.item():.2g}, BinnedRMSE: {binned_rmse_combo.item():.2g}, Flux: {flux_cost.item():.2g}, Eff: {model.efficiency.item():.6g}, Scale: {model.input_scale.item():.3g}, Loss: {loss_frac.item():.2g}')
            
            # Also save final storage and flux maps occasionally to watch learning progress
            with torch.no_grad():
                final_storage, total_flux, net_loss = model.compute_storage_and_flux()
                final_storage = final_storage.cpu().detach().numpy()
                total_flux = total_flux.cpu().detach().numpy()
                pd.DataFrame(final_storage).to_csv(f'{scenario_name}/EpochProgressMaps/Epoch{epoch}_FinalStorage.csv', header=False, index=False)
                pd.DataFrame(total_flux).to_csv(f'{scenario_name}/EpochProgressMaps/Epoch{epoch}_TotalFlux.csv', header=False, index=False)
            
            # Save history of model training progress
            with open(f'{scenario_name}/TrainingProgress.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Epoch', 'LearningRate', 'RMSE', 'RMSEbinned', 'Bias', 'FluxCost'])
                for epoch_i in range(epoch):
                    writer.writerow([epoch_i + 1, lr_list[epoch_i], rmse_list[epoch_i], rmse_binned_list[epoch_i],
                                     bias_list[epoch_i], flux_cost_list[epoch_i]])
        
        lr_list.append(current_lr)
        rmse_list.append(rmse_val.item())
        rmse_binned_list.append(binned_rmse_combo.item())
        bias_list.append(bias_val.item())
        flux_cost_list.append(flux_cost.item())
        
    # End of training, compute final results
    model.apply(clipper)
    final_storage, total_flux, net_loss = model.compute_storage_and_flux()
    
    # Parse weights into a more compact format
    outgoing_fracs = torch.zeros_like(final_storage)
    outgoing_fracs_up = torch.zeros_like(final_storage)
    outgoing_fracs_across = torch.zeros_like(final_storage)
    outgoing_fracs_down = torch.zeros_like(final_storage)
    for x, layer in enumerate(model.layers):
        outgoing_fracs[:,x] = layer.weight.data.sum(dim=1)[0:nY]
        for y in range(nY):
            if y > 0:
                outgoing_fracs_up[y,x] = layer.weight.data[y,y-1]
            outgoing_fracs_across[y,x] = layer.weight.data[y,y]
            if y < nY-1:
                outgoing_fracs_down[y,x] = layer.weight.data[y,y+1]
    
    final_storage = final_storage.cpu().detach().numpy()
    total_flux = total_flux.cpu().detach().numpy()
    net_loss = net_loss.cpu().detach().numpy()
    outgoing_fracs = outgoing_fracs.cpu().detach().numpy()
    outgoing_fracs_up = outgoing_fracs_up.cpu().detach().numpy()
    outgoing_fracs_across = outgoing_fracs_across.cpu().detach().numpy()
    outgoing_fracs_down = outgoing_fracs_down.cpu().detach().numpy()
    
    pd.DataFrame(final_storage).to_csv(f'{scenario_name}/FinalStorage.csv', header=False, index=False)
    pd.DataFrame(total_flux).to_csv(f'{scenario_name}/TotalFlux.csv', header=False, index=False)
    pd.DataFrame(net_loss).to_csv(f'{scenario_name}/NetLoss.csv', header=False, index=False)
    pd.DataFrame(outgoing_fracs).to_csv(f'{scenario_name}/TransferFracs.csv', header=False, index=False)
    pd.DataFrame(outgoing_fracs_up).to_csv(f'{scenario_name}/TransferFracs_Up.csv', header=False, index=False)
    pd.DataFrame(outgoing_fracs_across).to_csv(f'{scenario_name}/TransferFracs_Across.csv', header=False, index=False)
    pd.DataFrame(outgoing_fracs_down).to_csv(f'{scenario_name}/TransferFracs_Down.csv', header=False, index=False)
    
    input_scale = model.input_scale.cpu().detach().numpy()
    with open(f'{scenario_name}/InputScale.txt', 'w') as file:
        file.write(str(input_scale) + '\n')
    efficiency = model.efficiency.cpu().detach().numpy()
    with open(f'{scenario_name}/Efficiency.txt', 'w') as file:
        file.write(str(efficiency) + '\n')
    
    # Save history of model training progress
    with open(f'{scenario_name}/TrainingProgress.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Epoch', 'LearningRate', 'RMSE', 'RMSEbinned', 'Bias', 'FluxCost'])
        for epoch in range(n_epochs):
            writer.writerow([epoch + 1, lr_list[epoch], rmse_list[epoch], rmse_binned_list[epoch],
                             bias_list[epoch], flux_cost_list[epoch]])
    
#%%

if __name__ == '__main__':
  main()
