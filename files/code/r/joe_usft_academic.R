rm(list = ls())
library(tidyverse)
library(rvest)
library(lubridate)
library(tidylog)
library(ggplot2)
library(cowplot)
# Based on https://github.com/jlvoorheis/aea_webscrape/blob/master/aea_joe_compare_19_20.R

setwd("/Users/lukewatson/Documents/R_JOE/")

repec.t10 <- c("harvard","mit ","Massachusetts Institute of Technology","Princeton","Berkeley","University of Chicago",
               "Stanford","Columbia University","New York University","NYU","Yale","Brown")
repec.t25 <- c(repec.t10,"Boston University","University of California-San Diego","University of California San Diego","UCSD",
               "University of Michigan","University of Pennsylvania","UPenn","Dartmouth","Northwestern","Boston College",
               "University of Wisconsin-Madison", "University of Wisconsin Madison","University of California-Los Angeles",
               "University of California Los Angeles","UCLA","University of California Davis","University of California-Davis",
               "University of Southern California","Michigan State University","Cornell","Georgetown")

repec.t10 <- tolower(repec.t10)
repec.t25 <- tolower(repec.t25)

url.aea.19 <- "https://www.aeaweb.org/joe/listings?q=eNplj0FuwkAMRe_idZCAXXOASpW4QVVZZsYNQ41cJ7JcJ1QR4u41lLBhZz1_f_9_gV3xVnTw92pcJ-gvUBQptXJm6HUS6eCH599qGZ3J0gF6gA6cQ1IV-k_YwFcHR5ZFXtynuIXtevO2Wm9DXFytDEVJPl42qU7abEbj4e727-B05ozfVTKbLzCR5pKpMXoyOu3lmc84sTasKvOC5NEJIzcb7mMRvYK16RYbKX7LOAYlkZhHGh521-sf6V9eag,,"
aea_19 <- read_html(url.aea.19)
aea.19.txt <- data.frame(txt = read_lines(url.aea.19))
aea.19.txt$school.name <- str_detect(aea.19.txt$txt, "group-header-title")
aea.19.txt$post.date <- str_detect(aea.19.txt$txt, "listing-item-header-date-posted")
aea.19.txt <- aea.19.txt[aea.19.txt$school.name==T | aea.19.txt$post.date==T,]
aea.19.txt$txt <- qdapRegex::ex_between(aea.19.txt$txt, ">","<")
aea.19.txt$txt <- str_replace(aea.19.txt$txt, "Date Posted: ","")
aea.19.txt$school.name[aea.19.txt$school.name==T] <- aea.19.txt$txt[aea.19.txt$school.name==T]
aea.19.txt$school.name[aea.19.txt$school.name=="FALSE"] <- NA
aea.19.txt <- fill(aea.19.txt,school.name,.direction = "down")
aea.19.txt <- aea.19.txt[aea.19.txt$post.date==T,c("school.name","txt")]
aea.19.txt$school.name <- tolower(aea.19.txt$school.name)

aea.19.txt <- aea.19.txt %>%
  mutate(top.10 = str_detect(school.name,  paste(repec.t10, collapse="|")),
         top.25 = str_detect(school.name,  paste(repec.t25, collapse="|")))
aea.19.txt$date <-  as.Date(aea.19.txt$txt, format="%m/%d/%Y")
aea.19.txt$txt <- NULL


####

url.aea.20 <- "https://www.aeaweb.org/joe/listings?q=eNplj0Fqw0AMRe-itQuulz5AIZAbhFwilBnVnVTRGGmcYkLuXiWNV9mJ97--vq6wL96KTv5R7QzjFYpcIqVWLgwj9NDBD6-_1TI6k6XvgMGcw1EVxgO8w2cHXCcWGHUR6aC4L_fVoR_6t34Ic7UyFSXZvSipLtpsRePpkfaf4HThjF9VMptvMJHmkqkxejI6H4U3xTixNqwq64bk-RJGbzY8hhBvBWvLvTZS3JZ5DkpcIjHPND3jbrc_nQJdEw,,"
aea_20 <- read_html(url.aea.20)
aea.20.txt <- data.frame(txt = read_lines(url.aea.20))
aea.20.txt$school.name <- str_detect(aea.20.txt$txt, "group-header-title")
aea.20.txt$post.date <- str_detect(aea.20.txt$txt, "listing-item-header-date-posted")
aea.20.txt <- aea.20.txt[aea.20.txt$school.name==T | aea.20.txt$post.date==T,]
aea.20.txt$txt <- qdapRegex::ex_between(aea.20.txt$txt, ">","<")
aea.20.txt$txt <- str_replace(aea.20.txt$txt, "Date Posted: ","")
aea.20.txt$school.name[aea.20.txt$school.name==T] <- aea.20.txt$txt[aea.20.txt$school.name==T]
aea.20.txt$school.name[aea.20.txt$school.name=="FALSE"] <- NA
aea.20.txt <- fill(aea.20.txt,school.name,.direction = "down")
aea.20.txt <- aea.20.txt[aea.20.txt$post.date==T,c("school.name","txt")]
aea.20.txt$school.name <- tolower(aea.20.txt$school.name)

aea.20.txt <- aea.20.txt %>%
  mutate(top.10 = str_detect(school.name,  paste(repec.t10, collapse="|")),
         top.25 = str_detect(school.name,  paste(repec.t25, collapse="|")))
aea.20.txt$date <-  as.Date(aea.20.txt$txt, format="%m/%d/%Y")
aea.20.txt$txt <- NULL

####

aea.19.tot <- aea.19.txt %>%
  group_by(date,top.10) %>% 
  summarise(n_ads = n()) %>% 
  ungroup() %>% 
  group_by(top.10) %>% 
  mutate(total_ads = cumsum(n_ads), day_in_year = yday(date), year = as.factor(year(date)), month=month(date))

aea.20.tot <- aea.20.txt %>%
  group_by(date,top.10) %>% 
  summarise(n_ads = n()) %>% 
  ungroup() %>% 
  group_by(top.10) %>% 
  mutate(total_ads = cumsum(n_ads), day_in_year = yday(date), year = as.factor(year(date)), month=month(date))


all_aea <- bind_rows(aea.19.tot, aea.20.tot)

aea.19.txt.3 %>% 
  filter(month%in% c(8, 9, 10)) %>%
  ggplot() + 
  geom_line(aes(x=day_in_year, y=total_ads, group=top.10, colour=top.10))+
  scale_colour_brewer(palette="Set1")+theme_bw()

all_aea <- bind_rows(try19_df, try20_df)
all_aea$top_ten <- "top10"
all_aea$top_ten[all_aea$top.10==F] <- "nottop10"
all_aea$yr_top <- paste0(all_aea$year,"_",all_aea$top_ten)

# Color Blind Palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


gp <- all_aea %>% 
  filter(month%in% c(8, 9, 10)) %>%
  ggplot() +
  geom_line(aes(x=day_in_year, y=total_ads, colour=year, linetype=top_ten))+
  scale_colour_manual(values=cbp1) +
  scale_linetype_manual(values = c(rep("dashed", 1),rep("solid", 1))) +
  theme_bw() +
  labs(title="US FT Academic Cumulative Number of JOE Ads, 2019 vs. 2020",
       subtitle=paste0(" by REPEC Top 10; Date: ",Sys.Date()),
       x="Day in Year (Aug 1 = 213)",
       y="Cumulative Ads Posted on JOE",
       caption = "Based on John Voorheis") +
  theme(legend.position = "bottom") 
gp
save_plot("joe_us_ft_repec10.pdf", gp)
save_plot("joe_us_ft_repec10.png", gp)

####

aea.19.tot <- aea.19.txt %>%
  group_by(date,top.25) %>% 
  summarise(n_ads = n()) %>% 
  ungroup() %>% 
  group_by(top.25) %>% 
  mutate(total_ads = cumsum(n_ads), day_in_year = yday(date), year = as.factor(year(date)), month=month(date))

aea.20.tot <- aea.20.txt %>%
  group_by(date,top.25) %>% 
  summarise(n_ads = n()) %>% 
  ungroup() %>% 
  group_by(top.25) %>% 
  mutate(total_ads = cumsum(n_ads), day_in_year = yday(date), year = as.factor(year(date)), month=month(date))


all_aea <- bind_rows(aea.19.tot, aea.20.tot)
all_aea$top_25 <- "top25"
all_aea$top_25[all_aea$top.25==F] <- "nottop25"
all_aea$yr_top <- paste0(all_aea$year,"_",all_aea$top_25)

# Color Blind Palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


gp <- all_aea %>% 
  filter(month%in% c(8, 9, 10)) %>%
  ggplot() +
  geom_line(aes(x=day_in_year, y=total_ads, colour=year, linetype=top_25))+
  scale_colour_manual(values=cbp1) +
  scale_linetype_manual(values = c(rep("dashed", 1),rep("solid", 1))) +
  theme_bw() +
  labs(title="US FT Academic Cumulative Number of JOE Ads, 2019 vs. 2020",
       subtitle=paste0(" by REPEC Top 25; Date: ",Sys.Date()),
       x="Day in Year (Aug 1 = 213)",
       y="Cumulative Ads Posted on JOE",
       caption = "Based on John Voorheis") +
  theme(legend.position = "bottom") 
gp
save_plot("joe_us_ft_repec25.pdf", gp)
save_plot("joe_us_ft_repec25.png", gp)

