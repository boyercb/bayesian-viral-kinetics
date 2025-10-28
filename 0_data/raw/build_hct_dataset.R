library(haven)
library(tidyverse)
library(readxl)
library(pzfx)


# original hct ------------------------------------------------------------

sheets <- c(
  "qPCR nose",
  "qPCR throat",
  "FFA nose",
  "FFA throat"
)

hct <- map(sheets, function(sheet) {
  xl <- read_excel(
    path = "0_data/raw/COVHIC001_VL_qPCR_FFA_QCed_211109_SHARE.xlsx", 
    sheet = sheet,
    range = "A1:S39"
  )
  xl <- rename(xl, day = `Day post-inoculation/PID`)
  xl <- pivot_longer(xl, -day, names_to = "pid", values_to = "value")
  xl <- mutate(xl, name = str_replace(tolower(sheet), " ", "_"))
  xl
})

hct <- bind_rows(hct)
hct <- pivot_wider(hct, names_from = "name", values_from = "value")

hct <- arrange(hct, pid, day)

write_csv(hct, "0_data/hct_dat.csv")


# new  --------------------------------------------------------------------


# prism file
raw_pzfx <- "0_data/raw/Figure_1_PFU version 23.06.2022.pzfx"
  
# full symptom data
sx_dat <- read_excel("0_data/raw/symptoms HVO-CS-003 -Categorical SDC Data.xlsx")

# existing data
hct_dat <- read_csv("0_data/raw/hct_dat_old.csv")

# list all tables in Prism file
tabs <- pzfx_tables(raw_pzfx)

# get ids 
ids <- str_extract(tabs, "([0-9]{6})")
ids <- str_replace(ids, "667187", "667186")

# extract data from Prism files
dfs <- 
  map(1:length(tabs), 
      function(x) {
        df <- read_pzfx(raw_pzfx, table = x)
        if (x <= 18) {
          colnames(df) <- str_remove_all(colnames(df), "([0-9]{6})")
          colnames(df) <- tolower(colnames(df))
        } else if (x <= 36) {
          colnames(df) <- str_replace_all(colnames(df), "([0-9]{6})", "lfd")
          colnames(df)[1] <- "day"
        } else {
          colnames(df) <- str_replace_all(colnames(df), "([0-9]{6})", "symptom_score")
          colnames(df)[1] <- "day"
        }
        df$id <- ids[x]
        df
      })

# make nice names for symptoms
colnames(sx_dat) <- c(
  "quarantine",
  "cohort",
  "id",
  "challenge_time",
  "day",
  "day_desc",
  "timepoint",
  "timepoint_desc",
  "collection_method",
  "assessment_time",
  "runny_nose",
  "runny_nose_desc",
  "stuffy_nose",
  "stuffy_nose_desc",
  "sneezing",
  "sneezing_desc",
  "sore_throat",
  "sore_throat_desc",
  "hoarse_voice",
  "hoarse_voice_desc",
  "eye_soreness",
  "eye_soreness_desc",
  "earache",
  "earache_desc",
  "cough",
  "cough_desc",
  "chest_tightness",
  "chest_tightness_desc",
  "shortness of breath",
  "shortness of breath_desc",
  "wheeze",
  "wheeze_desc",
  "malaise",
  "malaise_desc",
  "headache",
  "headache_desc",
  "aches",
  "aches_desc",
  "fever",
  "fever_desc",
  "dizziness",
  "dizziness_desc",
  "rashes",
  "rashes_desc",
  "blisters",
  "blisters_desc",
  "diarrhoea",
  "diarrhoea_desc",
  "score"
)

# create intermediate data frames
pfus <- bind_rows(dfs[1:18])
lfd <- bind_rows(dfs[19:36])
sx_alt <- bind_rows(dfs[37:54])

# merge together
df <- left_join(pfus, lfd, by = c("id", "day"))
df <- left_join(df, sx_alt, by = c("id", "day"))

df$pid <- match(df$id, unique(ids))

# merge with existing data
hct_dat <- 
  left_join(hct_dat, df) |>
  relocate(day, id, pid) |>
  group_by(pid) |>
  mutate(id = first(id))

write_csv(hct_dat, "0_data/hct_dat.csv")
write_csv(sx_dat, "0_data/hct_sx_dat.csv")


ggplot(hct_dat, aes(x = time, y = rna, group = pid)) +
  facet_wrap(~factor(id, labels = unique(paste0("ID: ", pid)))) +
  geom_text(
    aes(
      x = time,
      y = 34,
      label = ifelse(lfd == 1, "X", "")
    ),
    color = "black",
    size = 2.5
  ) +
  geom_point(color = "#4ca5ff",  alpha = 0.75, shape = 16) +
  # geom_line(aes(x = time, y = rna_hat), color = "#4ca5ff") +
  geom_point(aes(x = time, y = pfu, group = id), color = "red", alpha = 0.5, shape = 16) +
  # geom_line(aes(x = time, y = pfu_hat), color = "red") +
  geom_point(aes(x = time, y = ffa, group = id), color = "green", alpha = 0.5, shape = 16) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal"
  ) +
  labs(
    x = "days from peak",
    y = "log count per ml"
  )
  