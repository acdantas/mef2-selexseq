# regex used for plotting

regex.l1 <- list( '01' = "^[CT]TA[TA]{4}TA[GA]$",              #.YTA(W4)TAR'
                  '02' = "^CTA[AT]{4}TAG",                     #.CTAW4TAG
                  '03' = "^TTA[AT]{4}TAA",                     #   TTAW4TAA
                  '04' = "(^CTA[TA]{4}TAA)|(^TTA[TA]{4}TAG)",  # CTAW4TAA
                  '05' = "^[CT]TA.*(G|C){1,4}.*TA[GA]$",       #YTA(noW4)TAR
                  '06' = "^(?![CT]TA)(.{3}[TA]{4}TA[GA])$",    #(noYTA)W4TAR 
                  '07' = "^[CT]TA[TA]{4}((?!TA[GA]).{3})$",    # .YTAW4(noTAR)
                  '08' = "^(?![CT]TA).{3}[TA]{4}((?!TA[GA]).{3})$"  # .(noYTA)W4(noTAR)
                 
)


re.l1.names <- list('01' = expression(YTA~W4~TAR),
                    '02' = expression(CTA~W4~TAG),
                    '03' = expression(TTA~W4~TAA),
                    '04' = expression(underline(CTA)~W4~underline(TAA)),
                    '05' = bquote("YTA"~underline("W4")~"TAR"),
                    '06' = bquote(underline("YTA")~"W4 TAR"),
                    '07' = bquote("YTA W4"~underline("TAR")),
                    '08' = bquote(underline("YTA")~"W4"~underline("TAR"))
                    
                    
)




regex.l2 <- list(CTAW4TAG = "CTA[AT][AT][AT][AT]TAG",
                 CTAW4TAA = "(CTA[AT][AT][AT][AT]TAA)|(TTA[AT][AT][AT][AT]TAG)",
                 TTAW4TAA = "TTA[AT][AT][AT][AT]TAA",                                            
                 'CTAW4(noTAR)'= "^CTA[TA]{4}(?!TA[GA]).{3}$",
                 '(noYTA)W4TAR'="^(?![CT]TA)(.{3}[TA]{4}TAG)$")

regex.l3 <- list('CTA(noW4)TAG' = "^CTA.*(G|C){1,4}.*TAG$",
                 'CTA(noW4)TAA' = "(^CTA.*(G|C){1,4}.*TAA$)|(^TTA.*(G|C){1,4}.*TAG$)",
                 'TTA(noW4)TAA' = "^TTA.*(G|C){1,4}.*TAA$",
                 'CTA(noW4/TAR)'= "^CTA.*(G|C){1,4}.*(?!TA[GA]).{3}$",
                 '(noYTA/W4)TAR'="^(?![CT]TA)(.{3}.*(G|C){1,4}.*TAG)$")


regex.l4 <- list(Primary = "CTA[ATCG][ATCG][ATCG][ATCG]TAG",
                                           Secondary = "(CTA[ATCG][ATCG][ATCG][ATCG]TA[AG])|([CT]TA[ATCG][ATCG][ATCG][ATCG]TAG)",
                                           Tertiary = "TTA[ATCG][ATCG][ATCG][ATCG]TAA") #, title="YTANNNNTAR")

regex.l5 <- list(Primary = "CTA[AT][AT][AT][AT]TAG",
                                           Secondary = "(CTA[AT][AT][AT][AT]TA[AG])|([CT]TA[AT][AT][AT][AT]TAG)",
                                           Tertiary = "TTA[AT][AT][AT][AT]TAA") #, title="YTAWWWWTAR")

regex.l6 <-  list(Primary = "CTA.*(G|C){1,4}.*TAG",
                                           Secondary = "(CTA.*(G|C){1,4}.*TA[AG])|([CT]TA.*(G|C){1,4}.*TAG)",
                                           Tertiary = "TTA.*(G|C){1,4}.*TAA") #, title="YTA(noW4tract)TAR")

regex.l7 <- list( 'YTA(W4)TAR'="^[CT]TA[TA]{4}TA[GA]$",
                  '4.CTAW4TAA' = "(^CTA[TA]{4}TAA)|(^TTA[TA]{4}TAG)",
                  '5.YTA(noW4)TAR'="^[CT]TA.*(G|C){1,4}.*TA[GA]$",
                  
                  '6.(noYTA)W4TAR'="^(?![CT]TA)(.{3}[TA]{4}TA[GA])$", 
                  '7.YTAW4(noTAR)'= "^[CT]TA[TA]{4}((?!TA[GA]).{3})$",  
                  '8.(noYTA)W4(noTAR)'= "^(?![CT]TA).{3}[TA]{4}((?!TA[GA]).{3})$"
                  
)

