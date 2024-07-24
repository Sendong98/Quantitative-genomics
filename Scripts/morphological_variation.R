library(readr)
birds <- read_delim("C:\\Users\\ggb24\\Downloads\\SequencedIndividualsBirdIDsExtra.txt")
birds$BirdID
View(birds)
plot(birds$Lifespan,birds$FROH)
abline(lm(birds$FROH~birds$Lifespan))
