# spoken-number-recognition
Program can distinguish between 3 spoken numbers from 1 speaker (can be trained for different speakers). After speech is recorded preprocessing is applied. Preprocessing contains silence and noise removal. After that LPC coefficents are calculated and used to extract features. Simple linear classificator is created which achieves high accuracy for 1 speaker. Train and test datasets used are also provided. Project also contains scripts for quantization of recorded signal and estimation of speaker's pitch frequency using two methods (autocorrelation method and parallel processing method).
