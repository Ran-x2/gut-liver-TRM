library(immunarch)
library(ggplot2)
immdata <- repLoad('F:/Tcell/immunarch_tissue/3')
div_chao <- repDiversity(immdata$data, "chao1")
p4 <- vis(div_chao)
p4

tc2 <- trackClonotypes(immdata$data, list("Donor AJD3280 LP", 10), .col = "aa+v")
p1 <- vis(tc2)
p1

immdata <- repLoad('F:/Tcell/immunarch_tissue/4')
div_chao <- repDiversity(immdata$data, "chao1")
p4 <- vis(div_chao)
p4

tc2 <- trackClonotypes(immdata$data, list("Donor AJG2309 LP", 10), .col = "aa+v")
p1 <- vis(tc2)
p1

immdata <- repLoad('F:/Tcell/immunarch_tissue/5')
div_chao <- repDiversity(immdata$data, "chao1")
p4 <- vis(div_chao)
p4

tc2 <- trackClonotypes(immdata$data, list("Donor AJKQ118 LP", 10), .col = "aa+v")
p1 <- vis(tc2)
p1


donors = c('AJD3280','AJG2309','AJKQ118')
for (i in c(3,4,5)){
  for (tissue in c('IEL','LP','L','PB')){
    name = paste0("Donor ",donors[i-2],' ',tissue)
    immdata <- repLoad(paste0('F:/Tcell/immunarch_tissue/', i))
    
    tc <- trackClonotypes(immdata$data, list(name, 10), .col = "aa+v")
    p <- vis(tc)
    ggsave(paste0(donors[i-2],"_",tissue,"_",".png"), plot = p, width = 6, height = 5, dpi = 300)
    
    
  }
}

donors = c('AJD3280','AJG2309','AJKQ118')
for (i in c(3,4,5)){
    immdata <- repLoad(paste0('F:/Tcell/immunarch_tissue/', i))
    imm_rare <- repClonality(immdata$data, .method = "rare")
    p <- vis(imm_rare) + vis(imm_rare, .by = "tissue", .meta = immdata$meta)
    ggsave(paste0(donors[i-2],"_rare.png"), plot = p, width = 12, height = 5, dpi = 300)
}
