library(readr)
library(readxl)
library(tidyverse)
library(lubridate)
birds <- read_delim("C:\\Users\\ggb24\\Downloads\\SequencedIndividualsBirdIDsExtra.txt")
birds$BirdID
View(birds)
vcf <- birds[birds$SeqID%in%filename$V1,]$BirdID
plot(birds$Lifespan,birds$FROH)
abline(lm(birds$FROH~birds$Lifespan))
birds_seq <- birds |> drop_na(Plate)
length(birds_seq$BirdID)
length(unique(birds_seq$BirdID))
body_mass <- read_xlsx("C:\\Users\\ggb24\\Downloads\\mass_fieldperiod.xlsx")
body_mass$day <-time_length(interval(body_mass$BirthDate,body_mass$OccasionDate),"day")
body_mass$month <-time_length(interval(body_mass$BirthDate,body_mass$OccasionDate),"month")

birth_mass1 <- body_mass[body_mass$month<=8,]
birth_mass1 <- birth_mass1 |> drop_na(BodyMass)
length(unique(birth_mass1$BirdID))*mean(unique(birth_mass1$BirdID)%in%birds_seq$BirdID)
length(unique(birth_mass1$BirdID))*mean(unique(birth_mass1$BirdID)%in%vcf)

birth_mass2 <- body_mass[body_mass$day<=7,]
birth_mass2 <- birth_mass2 |> drop_na(BodyMass)
length(unique(birth_mass2$BirdID))*mean(unique(birth_mass2$BirdID)%in%birds_seq$BirdID)
length(unique(birth_mass2$BirdID))*mean(unique(birth_mass2$BirdID)%in%vcf)
# sample size of birth mass: 100/99 (sequenced and birth mass data)
# sample size of juvenile mass: 1429/1369

# SNPs file sample id
filename <- read.table("C:\\Users\\ggb24\\Downloads\\filename.txt")

