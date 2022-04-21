library(tidyverse)
library(windowlickr)

# # Meta

meta = read_csv("data/metadata/sample-metadata.csv")

sidmel = meta %>%
    filter(SampleID%in%win$Samples,
           Species %in% c("Eucalyptus melliodora", "Eucalyptus sideroxylon"))


sidmel %>%
    select(SampleID, Species) %>%
    write_tsv("melsid.tsv")

# # extract a test window

bcf = "/data/kevin/work/workspaces/post-paneuc-big/data/variants-acanthophis/enhanced_filter/mpileup~bwa~Emelliodora_GCA_004368105.2~enhanced_filter.vcf.gz"

win = windowlickr:::bcf_getGT(bcf, region="CM023402.1:100000-110000")


pop_freq = NULL
for (i in seq_len(win$nSNP)) {
    f = tapply(win$GT[,i], sidmel$Species, function(x) {mean(x, na.rm=T)/2})
    pop_freq = cbind(pop_freq, f)
}

gt = win$GT

sid = meta %>%
    filter(Species=="Eucalyptus sideroxylon") %>%
    pull(SampleID)

mel = meta %>%
    filter(Species=="Eucalyptus melliodora") %>%
    pull(SampleID)

mic = meta %>%
    filter(Species=="Eucalyptus microcarpa") %>%
    pull(SampleID)
alb = meta %>%
    filter(Species=="Eucalyptus albens") %>%
    pull(SampleID)


sp1 = alb
sp2 = sid
dxy = NULL
for (l in seq_len(win$nSNP)) {
    d = 0
    n = 0
    for (i in match(sp1, win$Samples)) {
        for (j in match(sp2, win$Samples)) {
            I = gt[i, l]
            J = gt[j, l]
            if (is.na(I) || is.na(J)) next
            #print(paste(i, j, I, J))
            d = d + abs(I - J)
            n = n + 1
        }
    }
    dxy = c(dxy, l=d / n)
}

sum(dxy)/10000
mean(dxy)
hist(dxy)
plot(dxy)

