import pandas as pd
import numpy as np

electrode_df = pd.read_csv(#'Insert path to electrode coordinates csv file') #Ensure coordinate x column is titled "registered location.x'; same for y and z'


tissue_samples = [
    {'tissue_id': 1, 'x': , 'y': , 'z': },
    {'tissue_id': 2, 'x': , 'y': , 'z': },
    {'tissue_id': 3, 'x': , 'y': , 'z': },
    {'tissue_id': 4, 'x': , 'y': , 'z': },
    {'tissue_id': 5, 'x': , 'y': , 'z': },
    {'tissue_id': 6, 'x': , 'y': , 'z': },
    {'tissue_id': 7, 'x': , 'y': , 'z': },

    # Add more tissue sample coordinates here if needed
]


def euclidean_distance(x1, y1, z1, x2, y2, z2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


nearest_electrodes = []

for sample in tissue_samples:
    distances = electrode_df.apply(
        lambda row: euclidean_distance(row['registered location.x'], row['registered location.y'],
                                       row['registered location.z'], sample['x'], sample['y'], sample['z']), axis=1
    )

    sorted_indices = distances.nsmallest(4).index
    nearest_electrode_data = {
        'tissue_id': sample['tissue_id'],
        '1st_nearest_electrode': electrode_df.loc[sorted_indices[0], 'grid id'],
        '1st_distance': distances[sorted_indices[0]],
        '2nd_nearest_electrode': electrode_df.loc[sorted_indices[1], 'grid id'],
        '2nd_distance': distances[sorted_indices[1]],
        '3rd_nearest_electrode': electrode_df.loc[sorted_indices[2], 'grid id'],
        '3rd_distance': distances[sorted_indices[2]],
        '4th_nearest_electrode': electrode_df.loc[sorted_indices[3], 'grid id'],
        '4th_distance': distances[sorted_indices[3]]
    }
    nearest_electrodes.append(nearest_electrode_data)


nearest_electrode_df = pd.DataFrame(nearest_electrodes)


print(nearest_electrode_df.to_string())
