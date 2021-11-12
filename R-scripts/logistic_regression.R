
## Fit generalized linear model to frequency data ##

wd_glm = "" # <-- update here
setwd(wd_glm)

variant_props_fn = "" # <-- update here
variant_props = read.csv(variant_props_fn)
head(variant_props)

weeks = seq(1,nrow(variant_props), 1)

# here, var1_prop should be the weekly frequencies of a given variant
glm_var1 = glm(var1_prop ~ weeks, data = variant_props, family = binomial)
summary(glm_var1)

# visual inspection
plot(x = weeks, y = variant_props$var1_prop)
variant_props$predict.var1 = predict(glm_var1, type = 'response')
lines(variant_props$predict.var1)


