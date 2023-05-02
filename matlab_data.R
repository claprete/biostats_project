# Matlab data file
# 3/18/23
# Mac LaPrete

## Shigella strains over time ####
get_strains_matlab <- function(){
  times <- seq(2001,2020,length.out=20)
  # fraction of Shigella that is S.flexneri
  sflex <- c(0.72,0.85,0.72,0.72,0.78,0.79,0.80,0.74,0.80,0.78,
             0.79,0.85,0.64,0.77,0.59,0.73,0.54,0.66,0.75,0.56)
  # fraction of Shigella that is S.sonnei
  sson <- c(0.04,0.01,0.07,0.08,0.09,0.11,0.12,0.08,0.11,0.10,
            0.15,0.14,0.30,0.19,0.39,0.24,0.39,0.31,0.25,0.45)
  
  strains_matlab <- as.data.frame(cbind(times,sflex,sson))
  
  return(strains_matlab)
}

plot_strains_matlab <- function(strains_matlab){
  pdf("strains_matlab.pdf", width = 7, height = 5)
  plot(strains_matlab$times, strains_matlab$sflex,
       col="red", type = "l", ylim = c(0,1),
       main = "Matlab Site Strains Over Time",
       ylab = "Proportion of Shigella", xlab = "Years")
  lines(strains_matlab$times, strains_matlab$sson, col="blue")
  legend("topright",c("S. flexneri","S. sonnei"),col = c("red","blue"), lty = c(1,1))
  dev.off()
}


## Resistance over time ####
get_resistance_matlab <- function(){
  times <- c(2000,2005,2010,2015)
  
  cipro_all <- c(0.2, 14.5, 60.4, 83.6)/100
  cipro_sflex <- c(0, 17.3, 57.9, 82.8)/100
  cipro_sson <- c(0, 6.1, 86.5, 95.3)/100
  
  azith_all <- c(NA, NA, 27, 51.4)/100
  azith_sflex <- c(NA, NA, 13, 34)/100
  azith_sson <- c(NA, NA, 60, 91.3)/100
  
  mec_all <- c(2.7, 9.8, 12.3, 33.9)/100
  mec_sflex <- c(0.8, 11, 10.2, 35.7)/100
  mec_sson <- c(0, 2.7, 16.3, 31)/100
  
  cetr_all <- c(0, 0, 2.8, 12.2)/100
  cetr_sflex <- c(NA, 0, 2.8, 5.2)/100
  cetr_sson <- c(NA, 0, 3.4, 17.2)/100
  
  multidrug_all <- c(NA, 25.81, 29.23, 54.24)/100
  
  resistance_matlab <- as.data.frame(cbind(times,
                                          cipro_all, cipro_sflex, cipro_sson,
                                          azith_all, azith_sflex, azith_sson,
                                          mec_all, mec_sflex, mec_sson,
                                          cetr_all, cetr_sflex, cetr_sson,
                                          multidrug_all))
  
  return(resistance_matlab)
}

plot_resistance_matlab <- function(resistance_matlab){
  titles <- c("Ciprofloxacin All","Ciprofloxacin S. flexneri","Ciprofloxacin S. sonnei",
              "Azithromycin All","Azithromycin S. flexneri","Azithromycin S. sonnei",
              "Mecillinam All","Mecillinam S. flexneri","Mecillinam S. sonnei",
              "Ceftriaxone All","Ceftriaxone S. flexneri","Ceftriaxone S. sonnei",
              "Multidrug All")
  pdf("resistance_matlab.pdf", width = 12, height = 20)
  par(mfrow=c(5,3))
  for(i in 2:dim(resistance_matlab)[2]){
    plot(resistance_matlab$times, resistance_matlab[,2],
         col="red", type = "l", ylim = c(0,1),
         ylab = "Portion of Resistance", xlab = "Year",
         main = titles[i-1])
  }
  dev.off()
}






