# MMED project : Measles Transmission Dynamics DRC 2019

## Loading packages

library(tidyverse)
library(janitor)
library(lubridate)

## Read in data

cases <- read_csv("/Users/boazbaliejukia/Documents/PMF_Project/Measles_data_cleaning/MMEDGit/cases.csv") %>% clean_names()
coverage_who <- read_csv("/Users/boazbaliejukia/Documents/PMF_Project/Measles_data_cleaning/MMEDGit/Measles vaccination coverage DRC.csv") %>% clean_names()

## Data exploration

cases %>% head()
cases %>% glimpse()
coverage_who %>% head()

## Visualisation

#cases_18_19 <- cases %>% filter(year %in% c(2018, 2019)) 

cases <- cases %>% mutate(year_month <- ym(paste(year, month, sep = "-"))) %>% rename(
    "year_month" = `year_month <- ym(paste(year, month, sep = "-"))`
) #

cases %>% glimpse()

cases %>% filter(year %in% c(2019, 2020)) %>% 
    ggplot(aes(year_month, measles_suspect)) +
    geom_line(linewidth = 1, colour = BLUE) +
    scale_x_date(
        date_breaks = "2 months", 
        date_labels = "%b %Y"   # e.g. Jan 2019, Mar 2019
    ) +
    labs(
        title = "Democratic Republic of Congo (DRC)",
        subtitle = "Suspected measles cases during the 2019-2020 outbreak",
        x = "Time",
        y = "Suspected cases"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold",
                                  size = 24, color = "lightsalmon4"),
        plot.subtitle = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(t = 50)),
        
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.1),
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.1),
        text = element_text(size = 15),
        axis.ticks.length.y = unit(0, "mm"),
        axis.ticks.length.x = unit(2, "mm"),
        plot.margin = margin(0.05, 0.02, 0.05, 0.03, "npc")
    )

coverage_who$year <- as.character(coverage_who$year)
coverage_who %>% 
    filter(antigen == "MCV1" & 
               coverage_category == "OFFICIAL" & 
               year %in% c(2017, 2018, 2019, 2020)) %>% 
    ggplot(aes(year, coverage)) +
    geom_col(stat = "identity") +
    theme_bw() +
    labs(title = "Democratic Republic of Congo (DRC)",
         subtitle = "Routine Immunisation Coverage",
         x = "Year", y = "Coverage") +
    theme(
        plot.title = element_text(size = 24, face = "bold", colour = "lightsalmon4"),
        plot.subtitle = element_text(face = "bold", size = 20),
        legend.title = element_blank(),
        #panel.grid.major.y = element_line(color = "gray", linewidth = 0.1),
        text = element_text(size = 15),
        plot.margin = margin(0.05, 0.02, 0.05, 0.03, "npc")
    ) 

?geom_col()


