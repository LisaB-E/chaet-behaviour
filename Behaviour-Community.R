#' ---
#' title: “From behaviour to community change "
#' author: “Team KeithChaets"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: github_document
#' ---
#'

#' *Current model prepared by Gareth Last-name, under supervision by Dr Sally Keith.* 
#' *[LBE] Annotations = Lisa Bostrom Einarsson*
#' 
#' =======================================================================
#' # 1. Workspace preparation 
#' =======================================================================

rm(list=ls()) #removes all objects present in the workspace

#' Load packages
packages = c('reshape','gstat','ggplot2',
             'viridis','tictoc', 'progress')
load.pack = lapply(packages,require,char=T)
load.pack

#' =======================================================================
#' # 2. Preloading code
#' =======================================================================

#' ## IBM Parameters
ngenerations  = 10    # No. generations
replicates    = 10    # No. replicates
dim           = 100   # dimension of square habitat array
hab_dim       = dim^2 # total no. cells
nspecies      = 20    # No. species
nindiv        = 20    # No. individuals per species
tot_indiv     = nspecies*nindiv # Total individuals
eloss         = 30    # Time-step energy loss
fight_eloss   = 15    # Energy loss from aggression
exp_fight     = 20    # Roulette selection exponent
a             = 1
b             = 100
v             = 0.08  # Logistic mortality parameters
a_rep         = 0.4   # Asymptotic reproduction probability
offspring_pen = 2     # factor by which to scale offspring energy
diff_ag       = 0.05  # Differential aggression \in [0,0.5)

#' Potential one time-step movement. Can move to any adjecent cell (inc diagonal) 
step_moves = c(-1,dim-1,dim,dim+1,1,-dim+1,-dim,-dim-1)

#' Store species richness at each generation
rich = data.frame(gen=rep(1:ngenerations,replicates),
                  rep=rep(1:replicates,each=ngenerations),
                  rich=NA,
                  type=NA)
#'  Progress bar Windows
# only runs in windows!
# pb = winProgressBar(title='Generations',label='0%',min=0,max=100,initial=0,width=400)

#' Progress bar on Mac *use pb$tick() in for-loop to run*
pb = progress_bar$new(total=100, width=100, clear=F)

#' =======================================================================
#' #3. Functions
#' =======================================================================
#'
#' ## Aggress ------------------------------------------------------------

# Aggression relationships (transitive or intransitive) 
#' *LBE: I don't understand why intransitive (ie. non-linear) at smaller reps? Also, first time rps appears here - so throws an error*
if(reps<=replicates/2){
  intransitive = T 
}else{
  intransitive = F
}

#' Define aggression values
if(intransitive){
  aggression = matrix(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),ncol=nspecies)
  diag(aggression) = 0.5 # For interaction between individuals of the same species *OH, with intransitive vs transitive, does it mean intraspecific vs interspecific?*
}else{
  aggression = matrix(ncol=nspecies,nrow=nspecies,
                      sort(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),d=T))
}

#' ## Create habitat --------------------------------------------------------

## Simulate spatially autocorrelated habitat
#'  [LBE] The variogram is defined as the variance of the difference between field values at two locations vgm()= variogram. Exponential model. The experimental variogram is calculated by averaging onehalf the difference squared of the z-values over all pairs of observations with the specified separation distance and direction. In other words *the variogram describes textural differences between datasets, in terms of spatial locations, which common descriptive statistics and histograms do not incoporate*
hab_grid = expand.grid(1:dim, 1:dim) #A data frame containing one row for each combination of the supplied factors. The first factors vary fastest.
names(hab_grid) = c('x','y')
#' [LBE] look into whether these autocorrelated datasets are likley to be represenattive of coral reefs
hab_mod = gstat(formula=z~1,
                locations=~x+y, #ordinary & simple krigin
                dummy=T,
                beta=0,#Expected value of Gaussian field
                model=vgm(psill=0.025,model="Exp",range=30),# Range controls autocorrelation
                nmax=20)

hab_pred = predict(hab_mod, newdata=hab_grid, nsim=1) #predicts habitat value from model

#' Normalise to 0-1
hab_pred$sim1 = hab_pred$sim1 - min(hab_pred$sim1)
hab_pred$sim1 = hab_pred$sim1/max(hab_pred$sim1)

#' Make matrix
hab_vals = as.matrix(cast(hab_pred,
                          x~y,
                          value='sim1'))

#' ## Plot habitat values
hab.plot = ggplot(hab_pred) + theme_bw() +
  geom_raster(aes(x,y,fill=sim1)) +
  scale_x_continuous(expand=expand_scale(add=0)) +
  scale_y_continuous(expand=expand_scale(add=0)) +
  scale_fill_viridis(option='A',
                     name='Habitat values',
                     breaks = c(0,0.5,1),
                     guide=guide_colorbar(barheight = 0.5,
                                          barwidth = 10,
                                          ticks=F)) +
  theme(legend.position = 'top',
        axis.text=element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
#' 
#' ## Move --------------------------------------------------------------
#'
# delete me 
locs <- cbind(sample(hab_dim,tot_indiv),rep(1:nspecies,each=nindiv)) # Generate random starting positions for first run (or only to test code?)
habitat = matrix(0,ncol=dim,nrow=dim)
habitat[locs[,1]] = 1 
# delete me

 move <- function(locs, habitat){
   step_moves  <- c(-1,dim-1,dim,dim+1,1,-dim+1,-dim,-dim-1) #moving in all dimension (incl diagonal) [LBE]what about staying put?
   edge_corner <- unique(c(1:dim,                        # top row
                          which(1:hab_dim%%dim%in%0:1),  # sides
                          (hab_dim-dim+1):hab_dim))      # bottom
   ec_occ      <- edge_corner[which(habitat[edge_corner]==1)] # which edge corner cells are occupied
   ec_indiv    <- locs[locs[,1]%in%ec_occ,]                   # which individuals are on edge corner cells (V2&V1)
   offec_indiv <- locs[locs[,1]%in%setdiff(locs[,1],ec_indiv[,1]),] # individuals  NOT on edge corner cells (V1&v2)
   
   offec_indiv <- cbind(offec_indiv, sample(step_moves, dim(offec_indiv)[1], replace = TRUE)) # add move step (V3)

#Identify which edges or corners individuals on edges or corners are at...
#and define proposed moves
ec_indiv_new = NULL;
if(sum(ec_indiv[,1]==1)>0){                    #find if any inds are in top left corner
  top_left_ind = ec_indiv[ec_indiv[,1]==1,]    # if yes, 
  if(is.null(dim(top_left_ind))){
    top_left_ind = c(top_left_ind,sample(step_moves[3:5],1))
  }else{
    top_left_ind = cbind(top_left_ind,sample(step_moves[3:5],
                                             dim(top_left_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,top_left_ind)
}
if(sum(ec_indiv[,1] == hab_dim-dim+1)>0){
  top_right_ind = ec_indiv[ec_indiv[,1]==hab_dim-dim+1,]
  if(is.null(dim(top_right_ind))){
    top_right_ind = c(top_right_ind,sample(step_moves[5:7],1))
  }else{
    top_right_ind = cbind(top_right_ind,sample(step_moves[5:7],
                                               dim(top_right_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,top_right_ind)
}
if(sum(ec_indiv[,1] == dim)>0){
  bottom_left_ind = ec_indiv[ec_indiv[,1]==dim,]
  if(is.null(dim(bottom_left_ind))){
    bottom_left_ind = c(bottom_left_ind,sample(step_moves[1:3],1))
  }else{
    bottom_left_ind = cbind(bottom_left_ind,sample(step_moves[1:3],
                                                   dim(bottom_left_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,bottom_left_ind)
}
if(sum(ec_indiv[,1] == hab_dim)>0){
  bottom_right_ind = ec_indiv[ec_indiv[,1]==hab_dim,]
  if(is.null(dim(bottom_right_ind))){
    bottom_right_ind = c(bottom_right_ind,sample(step_moves[c(1,7:8)],1))
  }else{
    bottom_right_ind = cbind(bottom_right_ind,sample(step_moves[c(1,7:8)],
                                                     dim(bottom_right_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,bottom_right_ind)
}
if(sum(ec_indiv[,1]%in%(2:(dim-1)))>0){
  left_edge_ind = ec_indiv[ec_indiv[,1]%in%(2:(dim-1)),]
  if(is.null(dim(left_edge_ind))){
    left_edge_ind = c(left_edge_ind,sample(step_moves[1:5],1))
  }else{
    left_edge_ind = cbind(left_edge_ind,sample(step_moves[1:5],
                                               dim(left_edge_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,left_edge_ind)
}
if(sum(ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)))>0){
  right_edge_ind = ec_indiv[ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)),]
  if(is.null(dim(right_edge_ind))){
    right_edge_ind = c(right_edge_ind,sample(step_moves[c(1,5:8)],1))
  }else{
    right_edge_ind = cbind(right_edge_ind,sample(step_moves[c(1,5:8)],
                                                 dim(right_edge_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,right_edge_ind)
}
if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==1),c(1,hab_dim-dim+1)))>0){
  top_edge_ind = ec_indiv[ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==1),
                                                  c(1,hab_dim-dim+1)),]
  if(is.null(dim(top_edge_ind))){
    top_edge_ind = c(top_edge_ind,sample(step_moves[3:7],1))
  }else{
    top_edge_ind = cbind(top_edge_ind,sample(step_moves[3:7],
                                             dim(top_edge_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,top_edge_ind)
}
if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),c(dim,hab_dim)))>0){
  bottom_edge_ind = ec_indiv[ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),
                                                     c(dim,hab_dim)),]
  if(is.null(dim(bottom_edge_ind))){
    bottom_edge_ind = c(bottom_edge_ind,sample(step_moves[c(1:3,7:8)],1))
  }else{
    bottom_edge_ind = cbind(bottom_edge_ind,sample(step_moves[c(1:3,7:8)],
                                                   dim(bottom_edge_ind)[1],r=T))
  }
  ec_indiv_new = rbind(ec_indiv_new,bottom_edge_ind)
}

rownames(ec_indiv_new) = NULL

locs_new = rbind(offec_indiv,ec_indiv_new)
locs_new = locs_new[order(locs_new[,1]),]

# Habitat value of proposed new cell
locs_new = cbind(locs_new,hab_vals[locs_new[,1]+locs_new[,5]])

#Move probabilistically dependent upon the relative resource values of
#current and possible new cells
movers = locs_new[locs_new[,5]!=0,]
non_movers = locs_new[locs_new[,5]==0,]

move_vals = cbind(movers[,3],movers[,6])
rel_move_val = move_vals/rowSums(move_vals)
movers[runif(dim(movers)[1])<=rel_move_val[,1],5] = 0

locs_new = rbind(movers,non_movers)
locs_new = locs_new[order(locs_new[,1]),]

#Move
locs_new[,1] = locs_new[,1] + locs_new[,5]
locs_new[,3] = hab_vals[locs_new[,1]]
locs_new = locs_new[order(locs_new[,1]),]

#'
#' ## Die ---------------------------------------------------------------
#'
#'  Death occurs as Bernoulli trail with probability logistically dependent on current energy
mort = as.logical(rbinom(dim(locs)[1],1,a/(1+b*exp(-v*locs[,4]))))
locs = locs[mort,] #removes individuals from `locs` picked in `mort`

#' ## Reproduce ---------------------------------------------------------
#' ## Adjust energy -----------------------------------------------------
#' 
#' Loss energy,feed,die
locs[,4] = locs[,4]-eloss+locs[,3]*100; #loss of energy from one step, gain energy from habitat (feed)
locs[locs[,4]>100,4] = 100
#' 
#' =======================================================================
#' #4 Simulate
#' =======================================================================
#' 
#' #' Randomly place individuals of each species in habitat, record species ID,..
#' aggression, and habitat value
locs = cbind(sample(hab_dim,tot_indiv),rep(1:nspecies,each=nindiv)) # Randomly picks locations (V1-cell) for all individuals (V2) 
locs = cbind(locs,hab_vals[locs[,1]]) #V3-habitat value
#' add energetic value
locs = cbind(locs,rep(100,nindiv)) #v4 starting energetic value =100?
colnames(locs)=c("cell", "species", "hab_val", "energy")

#' Place individuals in habitat
habitat = matrix(0,ncol=dim,nrow=dim) #make 100 x 100 matrix
habitat[locs[,1]] = 1 # pulls out V1 (uniqe cell number) from locs, and places individual (=1) in the cell

#' Sort by location
locs = locs[order(locs[,1]),]


#' =======================================================================
#' #?. Not sure where this fits in
#' =======================================================================
#
#'  
  
  
  
  
#' 
#'
#' simulation start? *START here - think about breaking this up into functions instead? But do as branch in case you f it up...
 loop=0  
  for(generation in 1:ngenerations){
    loop = loop + 1
    #Progress bar update
    pb$tick()
   # time_update = proc.time() - start_time
   # info = sprintf('Generations run: %.2f%%. Time remaining: %.1f mins',
   #               100*loop/(ngenerations*replicates),
   #              time_update[3]/(loop*60)*(ngenerations*replicates-loop))
   # setWinProgressBar(pb,100*loop/(ngenerations*replicates),label=info)
    
    
    
    #Interacting species are those occupying the same cell
    loc_ind = locs_new[,1]
    int_loc = loc_ind[duplicated(loc_ind)]
    
    fighters = locs_new[locs_new[,1]%in%int_loc,]
    
    #Aggression is energetically expensive
    fighters[,4] = fighters[,4]-fight_eloss
    
    #Remove fighters from set of individuals
    locs_new = locs_new[locs_new[,1]%in%setdiff(locs_new[,1],fighters[,1]),]
    
    ###Roulette selection for winners
    all_winners = NULL
    
    for(i in unique(fighters[,1])){
      sp_in_cell = fighters[fighters[,1]==i,]
      sp_ids = sp_in_cell[,2]
      nosp = length(sp_ids)
      win_vect = rep(0,nosp)
      #Fights occur in random pairs, with the losers being removed until only one...
      #individual remains
      for(j in 1:(nosp-1)){
        fightj = sp_in_cell[sample(1:dim(sp_in_cell)[1],2),]
        
        aggress_vals = c(aggression[fightj[1,2],fightj[2,2]],
                         aggression[fightj[2,2],fightj[1,2]])
        
        energy_levels = c(fightj[1,4],fightj[2,4])
        
        status = aggress_vals*energy_levels
        #Status of individual is product of aggression and energy level. 
        #Highly aggressive but lower energy individuals can lose to less aggressive...
        #individuals who have high energy
        out_prob = cumsum(status^exp_fight)/sum(status^exp_fight)
        
        if(out_prob[1]>runif(1)){
          loser = fightj[2,]
        }else{
          loser = fightj[1,]
        }
        
        rowcount = 0
        rowsum = 0
        while(rowsum!=6){
          rowcount = rowcount + 1
          rowsum = sum(sp_in_cell[rowcount,] == loser)
        }
        sp_in_cell = sp_in_cell[-rowcount,]
      }
      
      initial = fighters[fighters[,1]==i,]
      
      rowcount = 0
      rowsum = 0
      while(rowsum!=6){
        rowcount = rowcount + 1
        rowsum = sum(initial[rowcount,] == sp_in_cell)
      }
      win_vect[rowcount] = 1
      all_winners = c(all_winners,win_vect)
    }
    
    fighters = cbind(fighters,all_winners)
    
    #Differentiate winners and losers
    winners = fighters[fighters[,7]==1,]
    losers = fighters[fighters[,7]==0,]
    
    ### Loss energy, consume, die
    #Losers only loss energy and die
    losers[,4] = losers[,4]-eloss
    losers = losers[losers[,4]>0,]
    mort = as.logical(rbinom(dim(losers)[1],1,a/(1+b*exp(-v*losers[,4]))))
    losers = losers[mort,]
    
    # Combines winner with those who didn't engage in aggression (i.e. were sole occupants of cell)
    locs_new = cbind(locs_new,all_winners=1)
    locs_new = rbind(locs_new,winners)
    locs_new = locs_new[order(locs_new[,1]),]
    
    #These individuals lose energy, consume, die and reproduce
    locs_new[,4] = locs_new[,4]-eloss+locs_new[,3]*100
    locs_new[locs_new[,4]>100,4] = 100
    locs_new = locs_new[locs_new[,4]>0,]
    mort = as.logical(rbinom(dim(locs_new)[1],1,a/(1+b*exp(-v*locs_new[,4]))))
    locs_new = locs_new[mort,]
    
    #Reproduction is Bernouilli trail dependent on current energy levels. 
    #Low-energy individuals (I sound like Trump!) are unlikely to reproduce.
    locs_new = cbind(locs_new,rbinom(dim(locs_new)[1],1,a_rep/(1+b*exp(-v*locs_new[,4]))))
    offspring = locs_new[locs_new[,8]==1,];
    #Set offspring energy to 50% of parents to account for lower chance ofsurvival
    offspring[,4] = offspring[,4]/offspring_pen;
    
    locs = rbind(locs_new[,1:4],offspring[,1:4],losers[,1:4])
    locs = locs[order(locs[,1]),]
    
    habitat = matrix(0,ncol=dim,nrow=dim);
    habitat[locs[,1]] = 1
    rich[rich$gen==generation & rich$rep==reps,3]=length(unique(locs[,2]))
    rich[rich$gen==generation & rich$rep==reps,4]=intransitive
    
  }
}
close(pb)

labs = c('Intransitive','Transitive')
count = 0
for(type in c(T,F)){
  count = count + 1
  subData = rich[rich$type==type,]
  
  frame = data.frame(gen=1:ngenerations,
                     mean=tapply(subData$rich,subData$gen,mean),
                     se=tapply(subData$rich,subData$gen,function(x) sd(x)/sqrt(length(x))),
                     type=labs[count])
  
  assign(labs[count],frame)
  
}


allRich = rbind(Intransitive,Transitive)

rich.plot = ggplot(allRich) + theme_bw() +
  geom_ribbon(aes(x=gen,ymin=mean-se,ymax=mean+se,fill=type),alpha=0.8) +
  scale_y_continuous(limits=c(0,nspecies)) +
  labs(y='Species richness',x='Generation') +
  scale_fill_manual(values=c("#35B779FF","#440154FF"),
                    name='Aggression type',
                    guide=guide_legend(keywidth = 1,
                                       keyheight=0.5)) +
  theme(legend.position='top',
        panel.grid=element_line(size=0.2),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7))

cairo_pdf('aggression.pdf',width=5,height=5.2)
print(rich.plot)
dev.off()
