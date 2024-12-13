# script to identify and list missing simulation runs with their seed numbers

library(stringr)

# running locally
directory_path = "/smdnas02/projects/jklus/JKSTProj/BlueHive_Sim_Results/SummaryLargePriorSS"

# running on server
# data_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/SummaryManuscript"
# directory_path = "/scratch/jklus/dpmSims/SummaryLargePriorSS"

# list file extensions
file_ext = list.dirs(directory_path)

# for each extension, make a new element of a data table with simulation attributes
summary_list = vector(mode = "list", length = length(file_ext)-1) 

for(ext_index in 2:length(file_ext)){ # skip first one, listing main directory
  # skip second one, Archive directory
  
  # initialize summary
  summary_table = data.frame(Path = NA, Model = NA, Scenario = NA, SM = NA, n_obs = NA, num_missing = NA) 
  
  # list files
  file_list = list.files(path = file_ext[ext_index])
 
  # parse out attributes and save
  short_file_ext = stringr::str_split(string = file_ext[ext_index], pattern = "/MODSUM")[[1]][2]
  parsed_dirname = unlist(stringr::str_extract_all(string = file_ext[ext_index], pattern = "_[:alnum:]+"))
  summary_table[1, "Path"] = file_ext[ext_index]
  summary_table[1, "Model"] = stringr::str_remove(string = parsed_dirname[3], pattern = "_")
  summary_table[1, "Scenario"] = stringr::str_remove(string = parsed_dirname[4], pattern = "_")
  summary_table[1, "n_obs"] = stringr::str_remove(string = parsed_dirname[5], pattern = "_n")
  summary_table[1, "SM"] = stringr::str_remove(string = parsed_dirname[6], pattern = "_")
  # temp_summary_table[1:length(file_list), "Tag"] = stringr::str_remove(string = parsed_dirname[length(parsed_dirname)], pattern = "_")
  
  
  
  if(length(file_list) > 0){
    
    # enumerate missing sims
    avail_sims = as.numeric(unlist(stringr::str_extract_all(string = file_list, pattern = "[:digit:]+")))
    missing_sims = which(!(1:100 %in% avail_sims))
    print(missing_sims)
    summary_table[1, "num_missing"] = length(missing_sims) 
    
    # check to catch any integer(0) returns from missingness check 
    if(length(missing_sims) == 0){
      missing_sims = 0
    } 

    
    summary_list[[ext_index-1]][[1]] = summary_table
    summary_list[[ext_index-1]][[2]] = missing_sims

  } else{
    
    # all sims missing
    summary_table[1, "num_missing"] = 100
    summary_list[[ext_index-1]][[1]] = summary_table
    summary_list[[2]] = 1:100
    
  }
  
}

summary_list[[1]]

# loop through and find which 