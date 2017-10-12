Instruction for using the Matlab code for reproducing results:
In order to create a signature:
Run runCreateSignature.m script in MATLAB environment.
It accepts two parameters: the folder where the data files reside, and a boolean parameter, where 0=African American training cohort and 1=European training cohort.
The output is a robust signature
To test the signature:
Run CompareSignature.m  in MATLAB environment.
It accepts two parameters: the folder where the data files reside, and a boolean parameter, where 0=African American validation cohort and 1=European validation cohort.
The output is the statistics (R2 and p-values)


