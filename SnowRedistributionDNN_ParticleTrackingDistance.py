# Copyright Eli Boardman 2024
# All rights reserved

import os
import torch
import pandas as pd
import csv

os.chdir('C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/')

device = torch.device('cuda')

#%%

class RedistributionNetwork(torch.nn.Module):
    def __init__(self, nX, nY,
                 outgoing_fracs_up, outgoing_fracs_across, outgoing_fracs_down,
                 efficiency):
        super(RedistributionNetwork, self).__init__()
        
        self.efficiency = torch.nn.Parameter(torch.tensor(efficiency).to(device))
        
        self.nX = nX
        self.nY = nY
        
        self.initial_flux =  torch.zeros(self.nY).to(device)
        
        self.layers = torch.nn.ModuleList()
        
        for x in range(nX):
            layer = torch.nn.Linear(nY, nY, bias=False)
            layer.weight.data.fill_(0)
            
            for y in range(nY):
                if y > 0:
                    layer.weight.data[y, y-1] = outgoing_fracs_up[y, x]
                layer.weight.data[y, y] = outgoing_fracs_across[y, x]
                if y < nY-1:
                    layer.weight.data[y, y+1] = outgoing_fracs_down[y, x]
            
            self.layers.append(layer)
    
    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
        return x
    
    def compute_storage(self, initial_storage, xStart):
        final_storage = torch.zeros_like(initial_storage)
        incoming_flux = self.initial_flux
        
        for x, layer in enumerate(self.layers):
            if x >= xStart:
                
                # Available mass from incoming flux and local initial storage
                total_avail = incoming_flux + initial_storage[:,x]
                
                # Remaining (final) mass storage
                final_storage[:,x] = total_avail * (1 - layer.weight.sum(dim=1))
                
                # Determine incoming transport to next layer
                incoming_flux = (total_avail @ layer.weight) * self.efficiency
                
        return final_storage

#%%

def load_data(file_path):
    data = pd.read_csv(file_path, header=None)
    return torch.tensor(data.values, dtype=torch.float64)

#%%

def main():
    
    scenario_name = 'Smoothed2x21_22.5deg_0.1loss'
    
    torch.set_default_dtype(torch.float64)
    
    # Load pre-trained flow fractions and initial pattern
    outgoing_fracs_up = load_data(f'{scenario_name}/TransferFracs_ExportOnly_Up.csv').to(device)
    outgoing_fracs_across = load_data(f'{scenario_name}/TransferFracs_ExportOnly_Across.csv').to(device)
    outgoing_fracs_down = load_data(f'{scenario_name}/TransferFracs_ExportOnly_Down.csv').to(device)
    
    net_export = load_data(f'{scenario_name}/NetExport.csv').to(device)
    
    with open(f'{scenario_name}/Efficiency.txt', 'r') as file:
        efficiency = float(file.read())
    
    nX = len(outgoing_fracs_across[0]) # Number of layers
    nY = len(outgoing_fracs_across)    # Number of neurons per layer
    
    model = RedistributionNetwork(nX, nY,
                                  outgoing_fracs_up,
                                  outgoing_fracs_across,
                                  outgoing_fracs_down,
                                  efficiency).to(device)
    
    # Track fate of one-hot inputs for each neuron
    onehot_input = torch.zeros_like(net_export).to(device)
    x_idxs = torch.arange(nX).unsqueeze(0).expand(nY, -1).to(device)
    y_idxs = torch.arange(nY).unsqueeze(1).expand(-1, nX).to(device)
    weighted_storage = torch.zeros_like(net_export).to(device)
    total_storage = torch.zeros_like(net_export).to(device)
    final_storage_stack = torch.zeros(nY, nX, nY, nX).to(device) # Caution - this can become very large!
    with torch.no_grad():
        for x in range(nX):
            for y in range(nY):
                if net_export[y, x] > 0:
                    onehot_input[y, x] = net_export[y, x]
                    
                    # Get map of hot cell's contribution to all other cells
                    final_storage = model.compute_storage(onehot_input, x)
                    final_storage_stack[y, x, :, :] = final_storage
                    total_storage = total_storage + final_storage
                    
                    # Weight relative to distance from hot cell
                    onehot_distance = torch.sqrt(torch.square(x_idxs - x) + torch.square(y_idxs - y))
                    weighted_storage = weighted_storage + onehot_distance * final_storage
                    
                    # For manual validation
                    #onehot_distance_data = onehot_distance.cpu().detach().numpy()
                    #pd.DataFrame(onehot_distance_data).to_csv(f'{scenario_name}/DistX{x}Y{y}.csv', header=False, index=False)
                    
                    onehot_input[y, x] = 0
            print(f'Done processing one-hot pattern for layer {x} of / {nX} total')
    
    storage_mask = total_storage > 0
    weighted_storage[storage_mask] = weighted_storage[storage_mask] / total_storage[storage_mask]
    
    weighted_storage_data = weighted_storage.cpu().detach().numpy()
    total_storage_data = total_storage.cpu().detach().numpy()
    pd.DataFrame(weighted_storage_data).to_csv(f'{scenario_name}/AvgTransportDist_ForNetImport.csv', header=False, index=False)
    pd.DataFrame(total_storage_data).to_csv(f'{scenario_name}/TotalImport_OneHotSum.csv', header=False, index=False)
    
    # Go through every neuron again and evaluate where its imported snow came from
    
    DistBinMax = 30
    DistBinSize = 1
    DistBins = list(range(0,DistBinMax,DistBinSize))
    nDist = len(DistBins)
    
    ContribThreshList = list([0.001,0.01,0.1])
    nContrib = len(ContribThreshList)
    
    transport_hist = torch.zeros(nY, nX, nDist).to(device)
    import_area = torch.ones(nY, nX, nContrib).to(device)
    for x in range(nX):
        for y in range(nY):
            if total_storage[y, x] > 0:
                
                # Distance from all other cells to current storage cell
                transport_distance = torch.sqrt(torch.square(x_idxs - x) + torch.square(y_idxs - y))
                
                for i, DistBin in enumerate(DistBins):
                    if DistBin == max(DistBins):
                        dist_mask = transport_distance >= DistBin
                    else:
                        dist_mask = (transport_distance >= DistBin) * (transport_distance < (DistBin + DistBinSize))
                    dist_mask = dist_mask * (net_export > 0)
                    
                    if dist_mask.any():
                        contrib_cells = torch.where(dist_mask)
                        transport_hist[y, x, i] = final_storage_stack[contrib_cells[0], contrib_cells[1], y, x].sum()
                
                # Contributing area
                for i, ContribThresh in enumerate(ContribThreshList):
                    contrib_mask = final_storage_stack[:, :, y, x] >= ContribThresh
                    import_area[y, x, i] += torch.count_nonzero(contrib_mask)
                    
        
        print(f'Done processing histogram for layer {x} of / {nX} total')
    
    transport_hist_data = transport_hist.cpu().detach().numpy()
    with open(f'{scenario_name}/TransportDistanceHistogram.csv', 'w', newline='') as csvfile:
        col_names = ['x','y','ExportTotal','ImportTotal','AvgDist']
        col_names = col_names + [f'ImportMinDist{DistBin}' for DistBin in DistBins]
        writer = csv.writer(csvfile)
        writer.writerow(col_names)
        for x in range(nX):
            for y in range(nY):
                row_data = transport_hist_data[y, x, :].tolist()
                writer.writerow([x, y,
                                net_export[y, x].item(),
                                total_storage_data[y, x],
                                weighted_storage_data[y, x],
                                *row_data])
            print(f'Done writing histogram for layer {x} of / {nX} total')
    
    import_area_data = import_area.cpu().detach().numpy()
    with open(f'{scenario_name}/ImportAreaMaps.csv', 'w', newline='') as csvfile:
        col_names = ['x','y']
        col_names = col_names + [f'ImportAreaMin{ContribThresh*1000:.0f}mm' for ContribThresh in ContribThreshList]
        writer = csv.writer(csvfile)
        writer.writerow(col_names)
        for x in range(nX):
            for y in range(nY):
                row_data = import_area_data[y, x, :].tolist()
                writer.writerow([x, y,
                                *row_data])
            print(f'Done writing area maps for layer {x} of / {nX} total')
    
    
#%%

if __name__ == '__main__':
  main()
