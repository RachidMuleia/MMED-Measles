# MMED project : Measles Transmission Dynamics DRC 2019

## Loading packages

library(tidyverse)
library(janitor)
library(lubridate)

## Read in data

data_dir <- if (dir.exists("Data")) "Data" else normalizePath("../Data")

cases <- read_csv(file.path(data_dir, "cases.csv"), show_col_types = FALSE) %>% clean_names()
coverage_who <- read_csv(file.path(data_dir, "Measles vaccination coverage DRC.csv"), show_col_types = FALSE) %>% clean_names()

## Data exploration

cases %>% head()
cases %>% glimpse()
coverage_who %>% head()
#cases %>% view()
## Visualisation

#cases_18_19 <- cases %>% filter(year %in% c(2018, 2019)) 

cases <- cases %>% mutate(year_month <- ym(paste(year, month, sep = "-"))) %>% rename(
    "year_month" = `year_month <- ym(paste(year, month, sep = "-"))`
) #

cases %>% glimpse()
cases_19 <- cases %>% filter(
    year %in% c(2019, 2020) #& month > 7
)

 


cases_19$month[cases_19$month == 1 & cases_19$year ==2020 ] <- 13
cases_19$month[cases_19$month == 2 & cases_19$year ==2020 ] <- 14
cases_19$month[cases_19$month == 3 & cases_19$year ==2020 ] <- 15
cases_19$month[cases_19$month == 4 & cases_19$year ==2020 ] <- 16
cases_19$month[cases_19$month == 5 & cases_19$year ==2020 ] <- 17
cases_19$month[cases_19$month == 6 & cases_19$year ==2020 ] <- 18
cases_19$month[cases_19$month == 7 & cases_19$year ==2020 ] <- 19
cases_19$month[cases_19$month == 8 & cases_19$year ==2020 ] <- 20
cases_19$month[cases_19$month == 9 & cases_19$year ==2020 ] <- 21
cases_19$month[cases_19$month == 10 & cases_19$year ==2020 ] <- 22
cases_19$month[cases_19$month == 11 & cases_19$year ==2020 ] <- 23
cases_19$month[cases_19$month == 12 & cases_19$year ==2020 ] <- 24

cases_19 %>% view()
cases_19_8 <- cases_19 %>% filter(
    month > 5
)
cases_19_8$month_2 <- 1:nrow(cases_19_8)

cases_19_8 %>% glimpse()

cases_plot <- cases %>% filter(year %in% c(2019, 2020)) %>% 
    ggplot(aes(year_month, measles_suspect)) +
    geom_line(linewidth = 1, colour = "blue") +
    scale_x_date(
        date_breaks = "2 months", 
        date_labels = "%b %Y"   # e.g. Jan 2019, Mar 2019
    ) +
    labs(
        #title = "Democratic Republic of Congo (DRC)",
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

cases_plot

coverage_who$year <- as.character(coverage_who$year)
coverage_who %>% 
    filter(antigen == "MCV1" & 
               coverage_category == "OFFICIAL" & 
               year %in% c(2017, 2018, 2019, 2020)) %>% 
    ggplot(aes(year, coverage)) +
    geom_col(stat = "identity") +
    theme_bw() +
    labs(#title = "Democratic Republic of Congo (DRC)",
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


