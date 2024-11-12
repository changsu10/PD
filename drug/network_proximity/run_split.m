Closet_Distance_ZScore_split('input_gene',0.5,800,900); 
Closet_Distance_ZScore_split('input_gene',0.5,901,1000); 
Closet_Distance_ZScore_split('input_gene',0.5,1001,1100); 
Closet_Distance_ZScore_split('input_gene',0.5,1101,1200); 
Closet_Distance_ZScore_split('input_gene',0.5,1201,1300); 
Closet_Distance_ZScore_split('input_gene',0.5,1301,1400); 
Closet_Distance_ZScore_split('input_gene',0.5,1401,1500); 
Closet_Distance_ZScore_split('input_gene',0.5,1501,1585); 


time matlab -nodesktop -nodisplay -r "run('run_split1.m'); exit;"&
time matlab -nodesktop -nodisplay -r "run('run_split2.m'); exit;"&
time matlab -nodesktop -nodisplay -r "run('run_split3.m'); exit;"&
time matlab -nodesktop -nodisplay -r "run('run_split4.m'); exit;"&
time matlab -nodesktop -nodisplay -r "run('run_split5.m'); exit;"&
time matlab -nodesktop -nodisplay -r "run('run_split6.m'); exit;"&
time matlab -nodesktop -nodisplay -r "run('run_split7.m'); exit;"&
time matlab -nodesktop -nodisplay -r "run('run_split8.m'); exit;"&