# Climate conceptual figure
library(dplyr)
library(lubridate)
library(ggplot2)

#convert temps from F to C
f2c <- function(x) {
  x<- (5/9) * (x-32)
  return(x)
}



clim<- read.csv("yourpath/Data Repo/climate_storrs_weather_station.csv")
head(clim, 20)
str(clim)

clim<-clim %>% mutate(Date = ymd(DATE), obs=f2c(TOBS), min = f2c(TMIN), max = f2c(TMAX))



# add in sample data
dates_for_graph <-read.csv("yourpath/Data Repo/sampling_dates.csv")

dates_for_graph$date<- ymd(dates_for_graph$date)

temps_on_dates<-data.frame(clim %>% filter(Date %in% dates_for_graph$date))



# What was the average temperature in the week preceding sampling for each sampling event?



temps_1wk_before_sample<-clim %>% filter(Date %in% dates_for_graph$date | Date %in% (dates_for_graph$date-1)| Date %in% (dates_for_graph$date-2)|
                                           Date %in% (dates_for_graph$date-3)| Date %in% (dates_for_graph$date-4)| Date %in% (dates_for_graph$date-5)|
                                           Date %in% (dates_for_graph$date-6)) %>%
  mutate(sample_period = c(rep(1,7),rep(3,7), rep(4,7),rep(6,7))) %>%
  group_by(sample_period) %>%
  summarize(meanT_week_before = mean(obs))


(sample_fig<-clim %>% filter(Date > "2019-09-14" & Date <"2020-04-13") %>% #show the 10 days before and after first and last samples taken
    ggplot(aes(x=Date, y = obs))+
    geom_ribbon(aes(ymin=min,ymax = max), fill = "lightgrey")+
    geom_line(color = "black", size = 2)+
    theme(panel.background = element_blank(), panel.grid = element_blank(),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          axis.text = element_text(size= 26))+
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y")+
    annotate(geom = "segment",
             x=as_date(dates_for_graph$date), xend=as_date(dates_for_graph$date),
             y=temps_on_dates$obs,yend = -14,
             arrow=arrow(),
             size=4, color=c("#925A44", "#155289","#489FA7","#335A30"))+
    geom_rect(aes(xmin = as_date(dates_for_graph$date)[1]-12.3, xmax = as_date(dates_for_graph$date)[1]+12.3,
                  ymin = temps_on_dates$obs[1]-3.4, ymax = temps_on_dates$obs[1]+3),
              fill = "white", color = "#925A44", size = 2)+
    geom_rect(aes(xmin = as_date(dates_for_graph$date)[2]-10, xmax = as_date(dates_for_graph$date)[2]+10,
                  ymin = temps_on_dates$obs[2]-3.4, ymax = temps_on_dates$obs[2]+3),
              fill = "white", color = "#155289", size = 2)+
    geom_rect(aes(xmin = as_date(dates_for_graph$date)[3]-10, xmax = as_date(dates_for_graph$date)[3]+10,
                  ymin = temps_on_dates$obs[3]-3.4, ymax = temps_on_dates$obs[3]+3),
              fill = "white", color = "#489FA7", size = 2)+
    geom_rect(aes(xmin = as_date(dates_for_graph$date)[4]-10, xmax = as_date(dates_for_graph$date)[4]+10,
                  ymin = temps_on_dates$obs[4]-3.4, ymax = temps_on_dates$obs[4]+3),
              fill = "white", color = "#335A30", size = 2)+
        annotate(geom="text",
             x=as_date(dates_for_graph$date),
             y=(temps_on_dates$obs),
             label=paste(round(temps_on_dates$obs,0),"°C" ,"\n(",
                         round(temps_1wk_before_sample$meanT_week_before,0),"°C)"),
             size= 9)+
    annotate(geom = "segment",
             x=c(rep(as_date("2019-10-30"),2),rep(as_date("2020-01-13"),2)),
             xend=c(rep(as_date("2019-11-10"),2),rep(as_date("2020-01-23"),2)),
             y=c(29,20,29,20),yend = c(29,20,29,20),
             size=16, color=c("#925A44", "#489FA7","#155289","#335A30"))+
    annotate(geom="text",
             x=c(as_date("2019-11-11"),as_date("2020-01-24"),as_date("2019-11-11"),as_date("2020-01-24")),
             y=c(29,29,20,20),
             label=c("Fall", "Winter","Emergence", "Spring"),
             hjust=0,
             size = 18)+
    labs(y="Temperature (°C)"))
sample_fig

# pdf("yourpath/Sampling_Figure.pdf",width = 15,height = 8.5)
# sample_fig
# dev.off()





