"beams" <-
function(eS, startvar){

  startvar=as.numeric(startvar)
  
  range=rbind(c(0,2*exp(startvar[1])), c(0,2*startvar[2]))
  
  tempsol=msecalc(eS, exp(startvar[1]), startvar[2], FALSE)[1];
  
  pp = rbind(c(exp(startvar[1]), startvar[2])
            ,c(exp(startvar[1]), startvar[2]))
         
  intv1 = exp(startvar[1])*(2/3);
  intv2 = startvar[2]*(2/3);
      
  bestsol = c(pp[1,1], pp[1,2], tempsol);
  secondbestsol = c(pp[2,1], pp[2,2], tempsol);
  
  iteration =0;
  
  while(1){
    
    for (c in 1:2){
      sp = rbind(c(pp[c,1],       pp[c,2]), 
                 c(pp[c,1]+intv1, pp[c,2]), 
                 c(pp[c,1]-intv1, pp[c,2]), 
                 c(pp[c,1] ,      pp[c,2]+intv2), 
                 c(pp[c,1],       pp[c,2]-intv2), 
                 c(pp[c,1]+intv1, pp[c,2]+intv2), 
                 c(pp[c,1]+intv1, pp[c,2]-intv2), 
                 c(pp[c,1]-intv1, pp[c,2]+intv2), 
                 c(pp[c,1]-intv1, pp[c,2]-intv2))
      sp=ifelse(sp<1,1,sp)
      for(a in 1:9){
        tempsol=msecalc(eS, sp[a,1], sp[a,2], FALSE)[1]
        if (tempsol < bestsol[3]){
          secondbestsol=bestsol
          bestsol=c(sp[a,1], sp[a,2],tempsol)
        }
      }
    }
  
    pp[1,] = c(bestsol[1],bestsol[2]);
    pp[2,] = c(secondbestsol[1],secondbestsol[2]);
    intv1 = intv1/2;
    intv2 = intv2/2;
    
    iteration = iteration + 1;
         
    #print(c(iteration, intv1, intv2, pp[1,], pp[2,], bestsol[3], secondbestsol[3]))
    if( (intv1 < 1 & intv2 < 1 ) & (secondbestsol[3]-bestsol[3]<1) ){
      break;
    }
  
  
  }

  return(c(bestsol[1:2]))
}

