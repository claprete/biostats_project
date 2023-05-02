# Dhaka data file
# 3/18/23
# Mac LaPrete

## Shigella strains over time ####
get_strains_dhaka <- function(){
  times <- seq(2001,2020,length.out=20)
  # fraction of Shigella that is S.flexneri
  sflex <- c(0.60,0.59,0.62,0.58,0.61,0.50,0.57,0.55,0.58,0.44,
             0.40,0.37,0.47,0.52,0.45,0.40,0.62,0.32,0.51,0.41)
  # fraction of Shigella that is S.sonnei
  sson <- c(0.10,0.10,0.10,0.20,0.10,0.20,0.26,0.10,0.11,0.30,
            0.40,0.47,0.40,0.41,0.48,0.51,0.25,0.52,0.44,0.48)
  
  strains_dhaka <- as.data.frame(cbind(times,sflex,sson))
  
  return(strains_dhaka)
}

plot_strains_dhaka <- function(strains_dhaka){
  pdf("strains_dhaka.pdf", width = 7, height = 5)
  plot(strains_dhaka$times, strains_dhaka$sflex,
       col="red", type = "l", ylim = c(0,1),
       main = "Dhaka Site Strains Over Time",
       ylab = "Proportion of Shigella", xlab = "Years")
  lines(strains_dhaka$times, strains_dhaka$sson, col="blue")
  legend("topright",c("S. flexneri","S. sonnei"),col = c("red","blue"), lty = c(1,1))
  dev.off()
}


## Resistance over time ####
get_resistance_dhaka <- function(){
  times <- c(2000,2005,2010,2015)
  
  cipro_all <- c(0.2, 18.4, 59.9, 72.3)/100
  cipro_sflex <- c(0, 23.1, 59.4, 72.8)/100
  cipro_sson <- c(0, 30.6, 83.8, 90)/100
  
  azith_all <- c(NA, 22.9, 24.1, 55.1)/100
  azith_sflex <- c(NA, 26.7, 20.3, 31.8)/100
  azith_sson <- c(NA, 26.7, 33.8, 81.8)/100
  
  mec_all <- c(0.5, 26.9, 16.7, 18.8)/100
  mec_sflex <- c(1, 43.3, 18.8, 12.9)/100
  mec_sson <- c(0, 0, 16.9, 23.3)/100
  
  cetr_all <- c(0, 2, 1.2, 8.3)/100
  cetr_sflex <- c(NA, 0, 0, 1.5)/100
  cetr_sson <- c(NA, 0, 2.9, 25.5)/100
  
  multidrug_all <- c(NA, NA, 27.03, 57.39)/100
  
  resistance_dhaka <- as.data.frame(cbind(times,
                                          cipro_all, cipro_sflex, cipro_sson,
                                          azith_all, azith_sflex, azith_sson,
                                          mec_all, mec_sflex, mec_sson,
                                          cetr_all, cetr_sflex, cetr_sson,
                                          multidrug_all))
  
  return(resistance_dhaka)
}

plot_resistance_dhaka <- function(resistance_dhaka){
  titles <- c("Ciprofloxacin All","Ciprofloxacin S. flexneri","Ciprofloxacin S. sonnei",
              "Azithromycin All","Azithromycin S. flexneri","Azithromycin S. sonnei",
              "Mecillinam All","Mecillinam S. flexneri","Mecillinam S. sonnei",
              "Ceftriaxone All","Ceftriaxone S. flexneri","Ceftriaxone S. sonnei",
              "Multidrug All")
  pdf("resistance_dhaka.pdf", width = 12, height = 20)
  par(mfrow=c(5,3))
  for(i in 2:dim(resistance_dhaka)[2]){
    plot(resistance_dhaka$times, resistance_dhaka[,2],
         col="darkgreen", type = "l", ylim = c(0,1),
         ylab = "Portion of Resistance", xlab = "Year",
         main = titles[i-1])
  }
  dev.off()
}
