###### **********************
## CLASS 9 - SOFTWARE TOOLS - SEQUENCE CLASSIFICATION AND RANDOM FOREST V1
##
## Luis Villasana
##
## 2024-10-04
##
##***************************


## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library(randomForest)
library(Biostrings)


# check your set-up

usethis::git_remotes()
# you should see an origin (in your github account) and upstream (my github account)

usethis::git_sitrep()

# _ Code your own solution -----

# CREATE A NEW BRANCH - do this after every pull request has been merged!

usethis::pr_init(branch = "formidable") # change formidable w/ more informative branch name!

usethis::git_sitrep()
# check that Git local (project) Current local branch is formidable
# and this is differet from the git_sitrep output from line38!

# now do your R work
# for this assignment, create a script file named
# <your-university-user-name-from-your-uog-email-address>.R
# for me, that would be cottenie.R
# this file will consist of your commented script file for the randomforest code
# put it in the correct folder!

# COMMIT your new file with the changes in the "git" tab next to Environmnet

# _ Submit pull request -----

usethis::pr_push(),

# in the launched browser window, click "Create pull request"
# in the description, add - "This is my random forest script file :-)"
# click "Create pull request"

# then you wait for me to review and accept - this is your "quiz" for this week

# after I have merged the pull request

usethis::pr_pull() # make sure your branch is up-to-date

usethis::pr_finish() # cleaning everything up, deleting local 

usethis::git_sitrep()


## load in data -------
df_te_t_spir = read_csv("../data/2026-10-26_t_spiralis.csv")
df_te_t_spir

counting_TE <- df_te_t_spir |> 
  count(matching_class) |> 
  unique() 

## _ second issue, from Livpreet's question, values in between parentheses ----
df_te_t_spir =  df_te_t_spir |> 
  mutate(qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
         # notice the ary_left instead of qry_left2 from line 73
         .after = qry_left) |> 
  mutate(in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
         .after = in_repeat_end)

## Two most commons TE ------
ls_non_te = c( "Simple_repeat",
               "Low_complexity", "Satellite",
               "rRNA", "snRNA", "ARTEFACT",
               "Unknown")
ls_non_te

# To observe the amount of repetitions of TE
df_te_t_spir |> 
  count(matching_class) |> 
  filter(matching_class %notin% ls_non_te) |> 
  arrange(desc(n))
##===========> The two most common TE matching_class groups are: LTR/Pao and PLE/Naiad

#To eliminate the Non_TE from the table
df_te_t_spir = df_te_t_spir |> 
  filter(matching_class %notin% ls_non_te)

# Checking it works
# none of the non-TE classes should remain
any(df_te_t_spir$matching_class %in% ls_non_te)   # should be FALSE
# remaining classes and how many rows each
df_te_t_spir |> 
  count(matching_class) |> 
  arrange(desc(n))

# Removing the other TE from the table, only leaving the top 2
ls_common_te = c("LTR/Pao", "PLE/Naiad")

df_te_t_spir = df_te_t_spir |> 
  filter(matching_class %in% ls_common_te)

# Checking if it needs cleaning
df_te_t_spir |> 
  mutate(n_bad = str_count(nucl_seq, "[^ACGTacgt]")) |> 
  summarise(rows_with_bad = sum(n_bad > 0),
            max_bad = max(n_bad))
#===========> It doesn´t need cleaning

# To obtain information
df_te_t_spir |> 
  mutate(len = nchar(nucl_seq)) |> 
  group_by(matching_class) |> 
  summarise(min = min(len), q1 = quantile(len, 0.25), median = median(len),
            q3 = quantile(len, 0.75), max = max(len))


#### 4 - CALCULATING SEQUENCE FEATURES ----

# Convert nucl_seq to a DNAStringSet (kept separate from the table)
df_te_t_spir = as.data.frame(df_te_t_spir)
seqs = DNAStringSet(str_to_upper(df_te_t_spir$nucl_seq))
class(seqs)

# Nucleotide counts (A, C, G, T only)
df_te_t_spir = cbind(df_te_t_spir,
                     as.data.frame(letterFrequency(seqs, letters = c("A", "C", "G", "T"))))

# Proportions of A, T and G (C is the remainder, so it adds no independent information)
total = df_te_t_spir$A + df_te_t_spir$C + df_te_t_spir$G + df_te_t_spir$T
df_te_t_spir$Aprop = df_te_t_spir$A / total
df_te_t_spir$Tprop = df_te_t_spir$T / total
df_te_t_spir$Gprop = df_te_t_spir$G / total

# Dinucleotide frequencies (k = 2, 16 features)
df_te_t_spir = cbind(df_te_t_spir,
                     as.data.frame(dinucleotideFrequency(seqs, as.prob = TRUE)))

# Trinucleotide frequencies (k = 3, 64 features)
df_te_t_spir = cbind(df_te_t_spir,
                     as.data.frame(trinucleotideFrequency(seqs, as.prob = TRUE)))



### 5 - TRAINING CLASSIFICATION MODEL: LTR/Pao vs PLE/Naiad ----

# Row id (stands in for processid) and a plain data frame
dfTE = df_te_t_spir |> 
  mutate(row_id = row_number()) |> 
  as_tibble()

table(dfTE$matching_class)

# Validation set: 20 per class, set aside until the end
set.seed(217)
dfValidation = dfTE |> 
  group_by(matching_class) |> 
  slice_sample(n = 20) |> 
  ungroup()

table(dfValidation$matching_class)

# Training set: rows not in validation, balanced between classes
n_train = dfTE |> 
  filter(row_id %notin% dfValidation$row_id) |> 
  count(matching_class) |> 
  pull(n) |> 
  min()
n_train   # size of the smaller class after removing validation

set.seed(13)
dfTraining = dfTE |> 
  filter(row_id %notin% dfValidation$row_id) |> 
  group_by(matching_class) |> 
  slice_sample(n = n_train) |> 
  ungroup()

table(dfTraining$matching_class)

# Feature sets, picked by name
f_comp  = c("Aprop", "Tprop", "Gprop")
f_di    = colnames(dinucleotideFrequency(DNAString("ACGT")))      
f_tri   = colnames(trinucleotideFrequency(DNAString("ACGT")))     

#nucleotide composition only
rf1 = randomForest(x = as.data.frame(dfTraining[, f_comp]),
                   y = factor(dfTraining$matching_class),
                   ntree = 500, importance = TRUE)
rf1

# Validation on unseen data
pred = predict(rf1, as.data.frame(dfValidation[, c(f_comp)]))


# Accuracy for the model
mean(pred == dfValidation$matching_class)