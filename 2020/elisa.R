setwd("E:/work/RProject/ELISA")
library(ggplot2)
library(data.table)
data <- read.csv("data.csv", stringsAsFactors = FALSE)
data <- read.csv("rawdata_new.csv", stringsAsFactors = FALSE)
data_no_outliner <- read.csv("data_no_outliner.csv", stringsAsFactors = FALSE)

data$treatment <- factor(data$treatment, levels = c("none", "control", "treat"))

windowsFonts(BL = windowsFont("黑体"))#設定字形
theme_set(theme_bw()+ theme(text = element_text(size=20,family = "BL"),legend.position="right"))#作圖設定
mycolors<-c("DDDDDD", "777777" ,"444444")

data <- data.table(data)
data_sum <- data[,.(BALP_mean = mean(BALP),BALP_sd = sd(BALP)/sqrt(10), BMP_2_mean = mean(BMP_2), BMP_2_sd = sd(BMP_2)/sqrt(10),
                    OCN_mean = mean(OCN),OCN_sd = sd(OCN)/sqrt(10), TGF_B1_mean = mean(TGF_B1),TGF_B1_sd = sd(TGF_B1)/sqrt(10))]


ggplot(data)+aes(treatment, TGF_B1, fill = week)+geom_boxplot()+labs(x = "分组", y = "TGF-B1(ng/ml)", legend())
ggplot(data)+aes(treatment, OCN, fill = week)+geom_boxplot()
ggplot(data)+aes(treatment, BALP, fill = week)+geom_boxplot()
ggplot(data)+aes(treatment, BMP_2, fill = week)+geom_boxplot()



ggsave(filename = "TGF-B1.png",
       ggplot(data)+aes(treatment, TGF_B1, fill = week)+geom_boxplot()+labs(x = "分组", y = "TGF-β1(ng/ml)")+
         scale_fill_grey(start = 0.9, end = 0.5,name="时间",labels=c("3周", "6周"))+
         scale_x_discrete(labels=c("空白组", "对照组", "实验组"))+scale_y_continuous(expand = c(0, 0),limits = c(0,19)),
       width = 20, height = 14, dpi = 300, units = "cm", device='png')
ggsave(filename = "OCN.png",
       ggplot(data)+aes(treatment, OCN, fill = week)+geom_boxplot()+labs(x = "分组", y = "OCN(ng/ml)")+
         scale_fill_grey(start = 0.9, end = 0.5,name="时间",labels=c("3周", "6周"))+
         scale_x_discrete(labels=c("空白组", "对照组", "实验组"))+scale_y_continuous(expand = c(0, 0),limits = c(0,1.2)),
       width = 20, height = 14, dpi = 300, units = "cm", device='png')
ggsave(filename = "BALP.png",
       ggplot(data)+aes(treatment, BALP, fill = week)+geom_boxplot()+labs(x = "分组", y = "BALP(ng/ml)")+
         scale_fill_grey(start = 0.9, end = 0.5,name="时间",labels=c("3周", "6周"))+
         scale_x_discrete(labels=c("空白组", "对照组", "实验组"))+scale_y_continuous(expand = c(0, 0),limits = c(0,26)),
       width = 20, height = 14, dpi = 300, units = "cm", device='png')
ggsave(filename = "BMP-2.png",
       ggplot(data)+aes(treatment, BMP_2, fill = week)+geom_boxplot()+labs(x = "分组", y = "BMP-2(ng/ml)")+
         scale_fill_grey(start = 0.9, end = 0.5,name="时间",labels=c("3周", "6周"))+
         scale_x_discrete(labels=c("空白组", "对照组", "实验组"))+scale_y_continuous(expand = c(0, 0),limits = c(0,1.8)),
       width = 20, height = 14, dpi = 300, units = "cm", device='png')





TGF_B1_aov <- aov(TGF_B1~treatment*week, data = data)
summary(TGF_B1_aov)
OCN_aov <- aov(OCN~treatment*week, data = data)
summary(OCN_aov)
BALP_aov <- aov(BALP~treatment*week, data = data)
summary(BALP_aov)
TukeyHSD(BALP_aov)
BMP_2_aov <- aov(BMP_2~treatment*week, data = data)
summary(BMP_2_aov)
TukeyHSD(BMP_2_aov)

sink("lm.csv")
print(summary(TGF_B1_aov))
sink()







ggplot(data_no_outliner)+aes(treatment, TGF_B1, fill = week)+geom_boxplot()
ggplot(data_no_outliner)+aes(treatment, OCN, fill = week)+geom_boxplot()
ggplot(data_no_outliner)+aes(treatment, BALP, fill = week)+geom_boxplot()
ggplot(data_no_outliner)+aes(treatment, BMP_2, fill = week)+geom_boxplot()


shapiro.test(data$BMP_2)
shapiro.test(data_no_outliner$TGF_B1)
shapiro.test(data_no_outliner$BALP)
shapiro.test(data_no_outliner$BMP_2)




TGF_B1_aov <- aov(TGF_B1~treatment*week, data = data_no_outliner)
summary(TGF_B1_aov)
OCN_aov <- aov(OCN~treatment*week, data = data_no_outliner)
summary(OCN_aov)
TukeyHSD(OCN_aov)
BALP_aov <- aov(BALP~treatment*week, data = data_no_outliner)
summary(BALP_aov)
TukeyHSD(BALP_aov)
BMP_2_aov <- aov(BMP_2~treatment*week, data = data_no_outliner)
summary(BMP_2_aov)
TukeyHSD(BMP_2_aov)


ggplot(data)+aes(treatment, TGF_B1.ng.ml., fill = week)+geom_boxplot()
ggplot(data)+aes(treatment, OCN.ng.ml., fill = week)+geom_boxplot()
ggplot(data)+aes(treatment, BALP.ng.ml., fill = week)+geom_boxplot()
ggplot(data)+aes(treatment, BMP_2.ng.ml., fill = week)+geom_boxplot()