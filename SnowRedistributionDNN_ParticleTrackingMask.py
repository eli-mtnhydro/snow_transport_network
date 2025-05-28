# Copyright Eli Boardman 2024
# All rights reserved

import os
import torch
import pandas as pd

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
    
    scenario_name = 'Smoothed2x21_22.5deg_0.01loss'
    mask_name = 'ContinentalGlacier'
    
    torch.set_default_dtype(torch.float64)
    
    # Load pre-trained flow fractions and initial pattern
    outgoing_fracs_up = load_data(f'{scenario_name}/TransferFracs_ExportOnly_Up.csv').to(device)
    outgoing_fracs_across = load_data(f'{scenario_name}/TransferFracs_ExportOnly_Across.csv').to(device)
    outgoing_fracs_down = load_data(f'{scenario_name}/TransferFracs_ExportOnly_Down.csv').to(device)
    
    net_export = load_data(f'{scenario_name}/NetExport.csv').to(device)
    mask_data = load_data(f'Mask_{mask_name}.csv').to(device)
    mask_inverse_data = torch.zeros_like(mask_data).to(device)
    mask_inverse_data[mask_data == 0] = 1
    
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
    mask_total = torch.zeros_like(net_export).to(device)
    mask_inverse_total = torch.zeros_like(net_export).to(device)
    with torch.no_grad():
        for x in range(nX):
            for y in range(nY):
                if net_export[y, x] > 0:
                    onehot_input[y, x] = net_export[y, x] * 1000
                    
                    # Get map of hot cell's contribution to all other cells
                    final_storage = model.compute_storage(onehot_input, x)
                    
                    # Calculate fraction of cell's total export that ends up in mask
                    mask_storage = mask_data * final_storage
                    mask_total[y, x] = mask_storage.sum() / 1000
                    
                    # Calculate fraction of cell's total export that ends up outside of mask
                    mask_inverse_storage = mask_inverse_data * final_storage
                    mask_inverse_total[y, x] = mask_inverse_storage.sum() / 1000
                    
                    onehot_input[y, x] = 0
            print(f'Done processing layer {x} of / {nX} total')
    
    mask_total_data = mask_total.cpu().detach().numpy()
    mask_inverse_total_data = mask_inverse_total.cpu().detach().numpy()
    pd.DataFrame(mask_total_data).to_csv(f'{scenario_name}/ExportContribToMask_{mask_name}.csv', header=False, index=False)
    pd.DataFrame(mask_inverse_total_data).to_csv(f'{scenario_name}/ExportExitFromMask_{mask_name}.csv', header=False, index=False)
    
#%%

if __name__ == "__main__":
  main()
