import pandas as pd
correct_cell = pd.read_csv("Baron_ActualCellType.csv") # actual labelled cell type for the unknown cells
cell_pred = pd.read_csv("cell_type_predictions_Baron.csv") # predictions

count = 0
count_unassigned = 0
for i in range(0, correct_cell.shape[0]):
    if correct_cell.iloc[i, 1] != cell_pred.iloc[i, 1]:
        count +=1
        if cell_pred.iloc[i, 1] == 'unassigned' and correct_cell.iloc[i, 1] != 'unassigned':
            count_unassigned +=1
        print("Cell Name: ", cell_pred.iloc[i,0])
        print("Actual: ", correct_cell.iloc[i, 1])
        print("Incorrect Prediction: ", cell_pred.iloc[i, 1], "\n")
        
print("Count differ: ", count)
print("Count of incorrect \"unassigned\" predictions: ", count_unassigned) 