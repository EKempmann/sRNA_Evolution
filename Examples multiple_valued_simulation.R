multi_value(mrna_name = "mrna_0_02_06.RData", 
            target_con = c(0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
            value_con = c(0.96, 0.97, 0.98, 0.985, 0.99), 
            target_gr = c(-13), 
            perc_of_mean_val = c(1.4),
            percentage_of_pos_val = c(0.25), 
            start_target_number_ = c(1),
            durchlaeufe = 100000,
            srna = "Ryhblong.fa",
            core = 8,
            wiederholungen_ = 2,
            time_steps_ = 10,
            store = 200,
            target_constrain = TRUE,
            rna_constrain = TRUE,
            accessibility_ = TRUE,
            keep_all_rnas_ = FALSE
            ) 

####### END OF COPY PASTE#######



# This command will produce 45 Output Files in the Current Working directory(9*5).
# In the Analyse graphics each point will represent core * wiederholungen_ Datapoints 

#ATTENTION ATTENTION!!!!!!!
#Name of the Outputfile can be confusing. Here is an example:
#valued_13_80_0_1.4_0.25_1_mrna_0_02_06.RData .RData
#
#        target_gr(-13 = 13), target_con(0.8 = 80), value_con(0.1=10), perc_of_mean_val, percentage_of_pos_val, start_target_number_, mRNA name
#valued_ 13         _         80          _         10      _          1.4      _       0.25       _           1           _          mrna_0_02_06.RData .RData
#(Yes the .RData occures 2times in the end because i was to lazy to crop it)
#If you change other Parameters as displayed in the name, 
# you may need to rename the data by hand or change the code responsible for the name (filename <- paste(...), easy to find)
#


#This Function  will prefor every possile combination of input setting.
# Parameters that will accept Vectors

#target_con: Target Value Constrain. 0.5 = 50%. Makes sure that the old Target Value wont decrease to much after each mutation checkpoint.
#value_con: Complete Value Constrain. 0.5 = 50%. Function compares Old sRNA value to new sRNA value and decide based on this parameter if the Mutations will be allowed.
# target_gr: Target threshold. Defines a Energie threshold. All mRNAs wich bounding Energy to the sRNA is below this threshold, will be seen as Targets.
#perc_of_mean_val: Defines the Value of the First sRNA Target based on the mean of all positive mRNA Values. 1.2 = 120% of the mean (WARNING: 10 = max value not 1000%)
# percentage_of_pos_val: Percentage of mRNAs with a positive value. 0.4 = 40% of the mRNAs will have a positive Value
# start_target_number_: Warning! this one can be irritating. Defines the number of Start Targets which will have a positive Value. All other Targets will have the Value 0.
# If the sRNA only has 4 Targets in the beginning, then 4 is the max you can type in here.


#Parameters that wont accept vecorts

#mrna_name: Data name of the mRNA file. Must be a RData file with the listed mRNA`s and their values as attributes. (See function value_for_mrna)
#durchlaeufe: Numbers of Mutations per mutation cycle
#srna: sRNA file (.fa or .txt but fasta format)
#core: specify the numbers of Cores that will be used. (Note: More Cores -> faster calculation but sometimes there occur errors if to many cores are used and you get 0byte Outputs)
#wiederholungen: times the Function will repeat itself. core * wiederholungen = number of trials
#time_steps: number of Mutation cycles that will be preformed. time_steps * durchlaeufe = Complete number of Mutations that will be inserted
#store: if 200 the function will check the sRNa and mRNA`s after 200 inserted Mutation and adopt or reject them, based on the Target Values.
#target_constrain / rna_constrain: FALSE equals target_con / value_con = 0. 
# accessibility: TRUE = intaRNA call will include the accessibility of the RNAs
# keep_all_rnas: FALSE = RNAs wont be safed in the output. WARNING If TRUE the Programm might crash. 
# This uses a LOT of RAM space (Because of my bad code) only use if you do small runs.
