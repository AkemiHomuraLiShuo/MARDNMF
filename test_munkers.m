A=[82,83,69,92;77,37,49,92;11,69,5,86;8,9,98,23];
[assignment,cost]=munkers(A)  
[assignedrows,dum]=find(assignment);
order=assignedrows'
r=munkers(A)  