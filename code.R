setwd("E:/work/Project/Dr.WWH_ELISA")
library(ggplot2)
library(broom)
library(data.table)

raw <- read.csv("data/2021_raw.csv")
windowsFonts(BL = windowsFont("黑体"))
theme_set(theme_bw() + theme(text = element_text(size = 20, family = "BL"),legend.position="bottom"))

#####label set####
raw$group <- factor(raw$group, levels = c("C", "L", "H"), labels = c("对照组", "低剂量", "高剂量"))
raw$week <- factor(raw$week, levels = c("3W", "6W"), labels = c("3周", "6周"))

####Mean and SE####
raw_sum <- data.table(raw)
all_sum <- raw_sum[,.(
    week = "All", group = "All",
    BALP_mean = mean(BALP),BALP_se = sd(BALP)/sqrt(nrow(raw_sum)), 
    BMP_2_mean = mean(BMP.2), BMP_2_se = sd(BMP.2)/sqrt(nrow(raw_sum)),
    IGF_1_mean = mean(IGF.1),IGF_1_se = sd(IGF.1)/sqrt(nrow(raw_sum))
    )]

#all_sum <- data.table(matrix(nrow = 0, ncol = 8))
#colnames(all_sum) <- c("week", "group", "BALP_mean", "BALP_se", "BMP_2_mean", "BMP_2_se", "IGF_1_mean", "IGF_1_se")

for(x in levels(raw_sum$week)){
    for(y in levels(raw_sum$group)){
        data <- raw_sum[week == x & group == y, ]
        data_sum <- data[, .(
            week = x, group = y,
            BALP_mean = mean(BALP), BALP_se = sd(BALP)/sqrt(nrow(data)), 
            BMP_2_mean = mean(BMP.2), BMP_2_se = sd(BMP.2)/sqrt(nrow(data)),
            IGF_1_mean = mean(IGF.1), IGF_1_se = sd(IGF.1)/sqrt(nrow(data))
        )]
        all_sum <- rbind(all_sum, data_sum)
    }
}

all_sum[, 3:8] <- round(all_sum[, 3:8], 2)
all_sum <- all_sum[, .(
    week = week, group = group,
    BALP = paste(BALP_mean, "±", BALP_se), BMP_2 = paste(BMP_2_mean, "±", BMP_2_se),  IGF_1 = paste(IGF_1_mean, "±", IGF_1_se)
)]
#fwrite(all_sum, "2021/all_sum.csv", row.names = FALSE)


####analysis function####
anl <- function(data, index){
    formula_week <- as.formula(paste0(index, " ~ week"))
    formula_group <- as.formula(paste0(index, " ~ group"))
    #different time
    result <- tidy(t.test(formula_week, data = data, var.equal = TRUE))
    #same group, different time
    result <- rbind(result, tidy(t.test(formula_week, data = data[data$group == "对照组",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_week, data = data[data$group == "低剂量",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_week, data = data[data$group == "高剂量",], var.equal = TRUE)))
    #different group
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "高剂量",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "低剂量",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "对照组",], var.equal = TRUE)))
    #same time, different group
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "高剂量" & data$week == "3周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "低剂量" & data$week == "3周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "对照组" & data$week == "3周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "高剂量" & data$week == "6周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "低剂量" & data$week == "6周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "对照组" & data$week == "6周",], var.equal = TRUE)))
    analysis <- c("time_all", "time_C", "time_L", "time_H", "group_all_C:L", "group_all_C:H", "group_all_L:H",
              "group_3W_C:L", "group_3W_C:H", "group_3W_L:H", "group_6W_C:L", "group_6W_C:H", "group_6W_L:H")
    result <- cbind(analysis, result)
    return(result)
}

####BALP####
ggsave(filename = "2021/BALP.png",
    ggplot(raw) + aes(x = group, y = BALP, fill = week) + geom_boxplot() +
        scale_y_continuous(expand = c(0,0), limits = c(0, 430)) + labs(x = "实验分组", y = "BALP(ng/ml)", fill = "饲养时间"),
    width = 20, height = 14, dpi = 600, units = "cm", device = "png")

result_BALP <- anl(raw, "BALP")
#write.csv(result_BALP, "2021/result_BALP.csv", row.names = FALSE)

####BMP-2####
ggsave(filename = "2021/BMP_2.png",
    ggplot(raw) + aes(x = group, y = BMP.2, fill = week) + geom_boxplot() +
        scale_y_continuous(expand = c(0,0), limits = c(0, 12.3)) + labs(x = "实验分组", y = "BMP-2(ng/ml)", fill = "饲养时间"),
    width = 20, height = 14, dpi = 600, units = "cm", device = "png")
    
result_BMP_2 <- anl(raw, "BMP.2")
#write.csv(result_BMP_2, "2021/result_BMP_2.csv", row.names = FALSE)

####IGF-1####
ggsave(filename = "2021/IGF_1_all.png",
    ggplot(raw) + aes(x = group, y = IGF.1, fill = week) + geom_boxplot() +
        scale_y_continuous(expand = c(0,0), limits = c(0, 29)) + labs(x = "实验分组", y = "IGF-1(ng/ml)", fill = "饲养时间"),#6W_C8: 27.93
    width = 20, height = 14, dpi = 600, units = "cm", device = "png")

ggsave(filename = "2021/IGF_1.png",    
    ggplot(raw) + aes(x = group, y = IGF.1, fill = week) + geom_boxplot() +
        scale_y_continuous(expand = c(0,0), limits = c(0, 10.5)) + labs(x = "实验分组", y = "IGF-1(ng/ml)", fill = "饲养时间"),
    width = 20, height = 14, dpi = 600, units = "cm", device = "png")    

result_IGF_1 <- anl(raw, "IGF.1")
#write.csv(result_IGF_1, "2021/result_IGF_1.csv", row.names = FALSE)