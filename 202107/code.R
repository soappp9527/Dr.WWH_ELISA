library(ggplot2)
library(broom)
library(data.table)

raw <- read.csv("202107/data.csv")
windowsFonts(BL = windowsFont("黑体"))
theme_set(theme_bw() + theme(text = element_text(size = 20, family = "BL"),legend.position="bottom"))

#####label set####
raw$group <- factor(raw$group, levels = c("空白组", "生理盐水", "金葡素"))
raw$week <- factor(raw$week, levels = c("3", "6"), labels = c("3周", "6周"))
raw[raw$IL.6.pg.ml. < 1, ]$IL.6.pg.ml. <- " < 1"

####Mean and SE####
raw <- data.table(raw)
all_sum <- raw[,.(
    week = "All", group = "All",
    BALP_mean = mean(BALP.ng.ml.),BALP_se = sd(BALP.ng.ml.)/sqrt(nrow(raw)), 
    BMP_2_mean = mean(BMP.2.ng.ml.), BMP_2_se = sd(BMP.2.ng.ml.)/sqrt(nrow(raw)),
    IL_6_mean = mean(IL.6.pg.ml.),IL_6_se = sd(IL.6.pg.ml.)/sqrt(nrow(raw))
)]

for(x in levels(raw$week)){
    for(y in levels(raw$group)){
        data <- raw[week == x & group == y, ]
        data_sum <- data[, .(
            week = x, group = y,
            BALP_mean = mean(BALP.ng.ml.),BALP_se = sd(BALP.ng.ml.)/sqrt(nrow(data)), 
            BMP_2_mean = mean(BMP.2.ng.ml.), BMP_2_se = sd(BMP.2.ng.ml.)/sqrt(nrow(data)),
            IL_6_mean = mean(IL.6.pg.ml.),IL_6_se = sd(IL.6.pg.ml.)/sqrt(nrow(data))
        )]
        all_sum <- rbind(all_sum, data_sum)
    }
}

all_sum[, 3:8] <- round(all_sum[, 3:8], 2)
all_sum <- all_sum[, .(
    week = week, group = group,
    BALP = paste(BALP_mean, "±", BALP_se), BMP_2 = paste(BMP_2_mean, "±", BMP_2_se),  IL_6 = paste(IL_6_mean, "±", IL_6_se)
)]
#fwrite(all_sum, "202107/all_sum.csv", row.names = FALSE)


####analysis function####
anl <- function(data, index){
    formula_week <- as.formula(paste0(index, " ~ week"))
    formula_group <- as.formula(paste0(index, " ~ group"))
    #different time
    result <- tidy(t.test(formula_week, data = data, var.equal = TRUE))
    #same group, different time
    result <- rbind(result, tidy(t.test(formula_week, data = data[data$group == "空白组",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_week, data = data[data$group == "生理盐水",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_week, data = data[data$group == "金葡素",], var.equal = TRUE)))
    #different group
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "金葡素",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "生理盐水",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "空白组",], var.equal = TRUE)))
    #same time, different group
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "金葡素" & data$week == "3周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "生理盐水" & data$week == "3周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "空白组" & data$week == "3周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "金葡素" & data$week == "6周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "生理盐水" & data$week == "6周",], var.equal = TRUE)))
    result <- rbind(result, tidy(t.test(formula_group, data = data[data$group != "空白组" & data$week == "6周",], var.equal = TRUE)))
    analysis <- c("time_all", "time_C", "time_S", "time_G", "group_all_C:S", "group_all_C:G", "group_all_S:G",
                  "group_3W_C:S", "group_3W_C:G", "group_3W_S:G", "group_6W_C:S", "group_6W_C:G", "group_6W_S:G")
    result <- cbind(analysis, result)
    return(result)
}


####BALP####
ggsave(filename = "202107/BALP.png",
       ggplot(raw) + aes(x = group, y = BALP.ng.ml., fill = week) + geom_boxplot() +
           scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + labs(x = "实验分组", y = "BALP(ng/ml)", fill = "饲养时间"),
       width = 20, height = 14, dpi = 600, units = "cm", device = "png")

result_BALP <- anl(raw, "BALP.ng.ml.")
#write.csv(result_BALP, "202107/result_BALP.csv", row.names = FALSE)

####BMP-2####
ggsave(filename = "202107/BMP_2.png",
       ggplot(raw) + aes(x = group, y = BMP.2.ng.ml., fill = week) + geom_boxplot() +
           scale_y_continuous(expand = c(0,0), limits = c(0, 7)) + labs(x = "实验分组", y = "BMP-2(ng/ml)", fill = "饲养时间"),
       width = 20, height = 14, dpi = 600, units = "cm", device = "png")

result_BMP_2 <- anl(raw, "BMP.2.ng.ml.")
#write.csv(result_BMP_2, "202107/result_BMP_2.csv", row.names = FALSE)

####IL-6####
ggsave(filename = "202107/IL_6.png",
       ggplot(raw) + aes(x = group, y = IL.6.pg.ml., fill = week) + 
           geom_hline(yintercept = 25, linetype = "dashed", color = "red", size = 0.6) +
           geom_dotplot(binwidth = 4.5, binaxis = "y", stackdir = "center", position = "dodge", method = "histodot") +
           scale_y_continuous(expand = c(0,0), limits = c(-5, 140), breaks = c(0, 25, 50, 75, 100, 125)) + 
           labs(x = "实验分组", y = "IL-6(pg/ml)", fill = "饲养时间"),
       width = 20, height = 14, dpi = 600, units = "cm", device = "png")

raw[, IL6_level := (ifelse(IL.6.pg.ml. < 25, "IL-6 ≤ 25", "IL-6 > 25"))]

IL6_group <- raw[, .N, by = c("week", "group", "IL6_level")]
IL6_group <- dcast.data.table(IL6_group, week + group ~ IL6_level, value.var = "N")
#write.csv(IL6_group, "202107/result_IL_6.csv", row.names = FALSE)