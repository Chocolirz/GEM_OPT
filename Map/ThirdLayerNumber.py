import math
import numpy as np

PlotThreeGEMsAligned = True
if PlotThreeGEMsAligned:
    # Your given drift results
    drift_results_aligned = np.array([
        4911, 
        2859, 2859, 2859, 2859, 2859, 2859, 
        902.7, 902.7, 902.7, 902.7, 902.7, 902.7, 
        511, 511, 511, 511, 511, 511, 
        82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 
        27.2, 27.2, 27.2, 27.2, 27.2, 27.2,
        4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 
        3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 
        0.7, 0.7, 0.7, 0.7, 0.7, 0.7
    ]) / 1e5

    # Lattice generation
    a = 0.0280
    R_max = 5 * a
    points = []
    max_rows = int(R_max / (math.sqrt(3)/2 * a)) + 2
    max_cols = int(R_max / a) + 2

    for row in range(-max_rows, max_rows + 1):
        y = row * (math.sqrt(3)/2) * a
        for col in range(-max_cols, max_cols + 1):
            x = col * a + (a/2 if row % 2 else 0)
            r = math.hypot(x, y)
            if r <= R_max + 1e-12:
                points.append((x, y, r))

    points.sort(key=lambda p: round(p[2], 12))
    positions = [(x, y) for x, y, _ in points]

    # First to second layer probability map
    prob_map = {pos: p for pos, p in zip(positions, drift_results_aligned)}

    # Convolution to get 1st → 3rd layer probability
    third_layer_map = {}
    for p1, prob1 in zip(positions, drift_results_aligned):
        for p2, prob2 in zip(positions, drift_results_aligned):
            x_new = p1[0] + p2[0]
            y_new = p1[1] + p2[1]
            third_layer_map[(x_new, y_new)] = third_layer_map.get((x_new, y_new), 0) + prob1 * prob2


double_values = np.round([x * 10100.2 for x in list(third_layer_map.values())[:19]],0)
photon_values = [int(x) for x in double_values]


# --------------------------------------------------------------------------------------------------------------------------------------------




PlotThreeGEMsMisaligned = True
if PlotThreeGEMsMisaligned:
    # Your given drift results
    drift_results_misaligned = np.array([
        4512.5, 3358, 1872.75, 4512.5, 1376, 3358, 1872.75, # r = 1/2, s3/2, s7/2, 1/2, 1.5, s3/2, s7/2
        769, 1872.75, 155, 1872.75, 155, 769, # r = s13/2, s7/2, s19/2, s7/2, s19/2, s13/2
        769, 122.25, 1376, 122.5, 769, 122.25, # s13/2, s21/2, 1.5, 2.5, s13/2, s21/2
        89.5, 57.25, 122.25, 18, 155, 16.75, 155, 16.75, 122.25, 18, 89.5, 57.25, 
        # 1.5s3, s31/2, s21/2, s37/2, s19/2, s39/2, s19/2, s39/2, s21/2, s37/2, 1.5s3, s31/2
        57.25, 7.5, 122.5, 5, 57.25, 7.5, # s31/2, s43/2, 2.5, 3.5, s31/2, s43/2
        5, 18, 0, 18, 0, 5, # 3.5, s37/2, _, s37/2, _, 3.5
        5, 0.75, 7.5, 0, 16.75, 0, 16.75, 0, 7.5, 0, 5, 0.75, 
        # 3.5, s57/2, s43/2, _, s39/2, _,           s39/2, _, s43/2, _, 3.5, s57/2
        0.75, 0, 5, 0, 0.75, 0 # s57/2, _, 3.5, _, s57/2, _
    ]) / 1e5

    # Lattice generation
    a = 0.0280
    R_max = 4 * a
    points = []
    max_rows = int(R_max / (math.sqrt(3)/2 * a)) + 2
    max_cols = int(R_max / a) + 2

    for row in range(-max_rows, max_rows + 1):
        y = row * (math.sqrt(3)/2) * a
        for col in range(-max_cols, max_cols + 1):
            x = col * a + (a/2 if row % 2 else 0)
            r = math.hypot(x, y)
            if r <= R_max + 1e-12:
                points.append((x, y, r))

    points.sort(key=lambda p: round(p[2], 12))
    positions = [(x + 0.5*a, y) for x, y, _ in points]
    
    # First to second layer probability map
    prob_map = {pos: p for pos, p in zip(positions, drift_results_misaligned)}

    # Convolution to get 1st → 3rd layer probability
    third_layer_map = {}
    for p1, prob1 in zip(positions, drift_results_misaligned):
        for p2, prob2 in zip(positions, drift_results_misaligned):
            x_new = p1[0] + p2[0]
            y_new = p1[1] + p2[1]
            third_layer_map[(x_new, y_new)] = third_layer_map.get((x_new, y_new), 0) + prob1 * prob2

double_values_inv = np.round([x * 10100.2 for x in list(third_layer_map.values())[:19]],0)
photon_values_inv = [int(x) for x in double_values_inv]