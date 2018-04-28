#' Bradley-Terry Weighted Likelihood Function with home parameter
#'
#' @param df Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' The row number is consistent with the team's index shown in dataframe. Column name must be "xn".
#' @param u The exponential decay rate.
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @return The estimated abilities.
#' @export

B_T_Weighted_home <- function(df,ability,u,i0){
  n1 <- nrow(df)
  n <- nrow(ability)-1
  stop <- 0
  j <- 1
  ability[,1] <- 0
  ability0 <- matrix(0,nrow = n,ncol = 2)

  for(i in 1:n1){
    ability0[df[i,1],1] <- ability0[df[i,1],1] + df[i,3]
    ability0[df[i,2],1] <- ability0[df[i,2],1] + df[i,4]
    ability0[df[i,1],2] <- ability0[df[i,1],2] + 1
    ability0[df[i,2],2] <- ability0[df[i,2],2] + 1
  }

  k1 <- c()
  k0 <- c()
  for(i in 1:n){
    if(ability0[i,1]==ability0[i,2]){
      k1 <- c(k1,i)
    } else if(ability0[i,1]==0){
      k0 <- c(k0,i)
    }
  }

  while(stop==0){
    F <- matrix(0,nrow = (n+1),ncol = 1)
    F1 <- matrix(0,nrow = (n+1),ncol = 1)
    ##s1 <- 0
    for(i in 1:n1){
      a1 <- df[i,1]
      a2 <- df[i,2]
      C <- exp(-u*df[i,5])
      x <- ability[a1,"xn"]
      y <- ability[a2,"xn"]
      p <- ability["at.home","xn"]+x-y
      q <- exp(p)
      A <- -(df[i,3]*(1/(q+1))+df[i,4]*(-q/(q+1)))*C
      B <- -(df[i,3]*(-q/(q+1)^2)+df[i,4]*(-q/(q+1)^2))*C
      F[a1,1] <- F[a1,1] + A
      F[a2,1] <- F[a2,1] - A
      F[n+1,1] <- F[n+1,1] + A
      F1[a1,1] <- F1[a1,1] + B
      F1[a2,1] <- F1[a2,1] + B
      F1[n+1,1] <- F1[n+1,1] + B
      ##s1 <- s1 - C*df[i,3]*(p-log(q+1))+df[i,4]*(-log(q+1))
    }

    F[,1] <- F[,1]-F[i0,1]
    F1[,1] <- 1/F1[,1]
    if(!identical(k0,NULL)){
      F[k0,1] <- 0
      if(ability0[i0,1]==0){
        ability[k0,1] <- 0
      } else{
        ability[k0,1] <- -1
      }
    }
    if(!identical(k1,NULL)){
      F[k1,1] <- 0
      ability[k1,1] <- 4
    }
    ability[,"xn"] <- ability[,"xn"]-F1*F/(1/(j-0.9)+1)

    ##cat(k,ability[21,1],'\n')
    k <- sum(abs(F[,1]))

    j <- j + 1
    if((k<0.00001)||(j>10000)){
      stop <- 1
    }
    cat(k,'\n')
  }
  return(ability)
}

#' Bradley-Terry Weighted Likelihood Function with home parameter for Likelihood Ratio Test of
#' abilities time constancy
#'
#' @param df Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter  whose rowname must be "at.home".
#' The row number is consistent with the team's index shown in dataframe. Column name must be "xn".
#' @param uf Team index whose time constancy will be tested.
#' @param u The exponential decay rate.
#' @param i0 is a teams index whose ability will be fixed as 0 (usually the team loss most).
#' @param date the date were selected team will have different abilities before and after
#' (usually chosen to be the mid time when abilities are expected to change).
#' @return The estimated abilities.
#' @export

B_T_Weighted_fp <- function(df,ability,uf,u,i0,date){
  n1 <- nrow(df)
  n <- nrow(ability)-1
  stop <- 0
  j <- 1
  n2 <- length(date)-1
  ability0 <- matrix(0,nrow = (n+1),ncol = n2)
  rownames(ability0) <- rownames(ability)
  colnames(ability0) <- seq(1,n2,1)

  while(stop==0){
    F <- matrix(0,nrow = (n+1),ncol = 1)
    F0 <- matrix(0,nrow = n2,ncol = 1)
    F1 <- matrix(0,nrow = (n+1),ncol = 1)
    F01 <- matrix(0,nrow = n2,ncol = 1)
    ##s1 <- 0
    for(i in 1:n1){
      a1 <- df[i,1]
      a2 <- df[i,2]
      C <- exp(-u*df[i,5])
      g <- 0

      if(a1==uf){
        stop1 <- 0
        m <- 2
        while(stop1==0){
          if(m==(n2+1)){
            x <- ability0[uf,n2]
            g <- m-1
            stop1 <- 1
          } else if(i<date[m]){
            x <- ability0[uf,(m-1)]
            g <- m-1
            stop1 <- 1
          } else{
          }
          m <- m + 1
        }
      } else{
        x <- ability0[a1,1]
      }

      if(a2==uf){
        stop1 <- 0
        m <- 2
        while(stop1==0){
          if(m==(n2+1)){
            y <- ability0[uf,n2]
            g <- m-1
            stop1 <- 1
          } else if(i<date[m]){
            y <- ability0[uf,(m-1)]
            g <- m-1
            stop1 <- 1
          } else{
          }
          m <- m + 1
        }
      } else{
        y <- ability0[a2,1]
      }
      p <- ability0["at.home",1]+x-y
      q <- exp(p)
      A <- -(df[i,3]*(1/(q+1))+df[i,4]*(-q/(q+1)))*C
      B <- -(df[i,3]*(-q/(q+1)^2)+df[i,4]*(-q/(q+1)^2))*C
      F[a1,1] <- F[a1,1] + A
      F[a2,1] <- F[a2,1] - A
      F[n+1,1] <- F[n+1,1] + A
      F1[a1,1] <- F1[a1,1] + B
      F1[a2,1] <- F1[a2,1] + B
      F1[n+1,1] <- F1[n+1,1] + B
      if(a1==uf){
        F0[g,1] <- F0[g,1] + A
        F01[g,1] <- F01[g,1] + B
      }
      if(a2==uf){
        F0[g,1] <- F0[g,1] - A
        F01[g,1] <- F01[g,1] + B
      }
      ##s1 <- s1 - C*df[i,3]*(p-log(q+1))+df[i,4]*(-log(q+1))
    }

    F[,1] <- F[,1]-F[i0,1]
    F1[,1] <- 1/F1[,1]
    F01[,1] <- 1/F01[,1]

    ability0[-uf,] <- ability0[-uf,]-F1[-uf,]*F[-uf,]/(1/(j-0.9)+1)
    ability0[uf,] <- ability0[uf,]-t(F01*F0/(1/(j-0.9)+1))

    ##cat(k,ability[21,1],'\n')
    k <- sum(abs(F[-uf,1]))+sum(abs(F0[,1]))
    ##print(k)
    j <- j + 1
    if((k<0.00001)||(j>10000)){
      stop <- 1
    }
  }
  return(ability0)
}

#' Bradley-Terry Weighted Likelihood Function with fixed home parameter
#'
#' @param df Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' The row number is consistent with the team's index shown in dataframe. Column name must be "xn".
#' @param home
#' @param u The exponential decay rate.
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @return The estimated abilities.
#' @export

B_T_Weighted_fhome <- function(df,ability,home,u,i0){
  n1 <- nrow(df)
  n <- nrow(ability)-1
  stop <- 0
  j <- 1

  ability <- as.matrix(ability[1:n,])
  colnames(ability) <- c("xn")

  while(stop==0){
    F <- matrix(0,nrow = n,ncol = 1)
    F1 <- matrix(0,nrow = n,ncol = 1)
    ##s1 <- 0
    for(i in 1:n1){
      a1 <- df[i,1]
      a2 <- df[i,2]
      C <- exp(-u*df[i,5])
      x <- ability[a1,"xn"]
      y <- ability[a2,"xn"]
      p <- home+x-y
      q <- exp(p)
      A <- -(df[i,3]*(1/(q+1))+df[i,4]*(-q/(q+1)))*C
      B <- -(df[i,3]*(-q/(q+1)^2)+df[i,4]*(-q/(q+1)^2))*C
      F[a1,1] <- F[a1,1] + A
      F[a2,1] <- F[a2,1] - A
      F1[a1,1] <- F1[a1,1] + B
      F1[a2,1] <- F1[a2,1] + B
      ##s1 <- s1 - C*df[i,3]*(p-log(q+1))+df[i,4]*(-log(q+1))
    }

    F[,1] <- F[,1]-F[i0,1]
    F1[,1] <- 1/F1[,1]

    ability[,"xn"] <- ability[,"xn"]-F1*F/(1/(j-0.9)+1)

    ##cat(k,ability[21,1],'\n')
    k <- sum(abs(F[,1]))

    j <- j + 1
    if((k<0.00001)||(j>10000)){
      stop <- 1
    }
  }
  ability <- rbind(ability,matrix(home))
  rownames(ability)[n+1] <- c("at.home")
  return(ability)
}

#' Bradley-Terry Weighted Likelihood Function with home parameter and return Confidence Interval
#' of estimated parameters
#'
#' @param df Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' The row number is consistent with the team's index shown in dataframe. Column name must be "xn".
#' @param u The exponential decay rate.
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @return The estimated abilities.
#' @export

B_T_Weighted_home_CI <- function(df,ability,u,i0){
  n1 <- nrow(df)
  n <- nrow(ability)-1
  stop <- 0
  j <- 1
  ability[,1] <- 0

  while(stop==0){
    F <- matrix(0,nrow = (n+1),ncol = 1)
    F1 <- matrix(0,nrow = (n+1),ncol = 1)
    ##s1 <- 0
    for(i in 1:n1){
      a1 <- df[i,1]
      a2 <- df[i,2]
      C <- exp(-u*df[i,5])
      x <- ability[a1,"xn"]
      y <- ability[a2,"xn"]
      p <- ability["at.home","xn"]+x-y
      q <- exp(p)
      A <- -(df[i,3]*(1/(q+1))+df[i,4]*(-q/(q+1)))*C
      B <- -(df[i,3]*(-q/(q+1)^2)+df[i,4]*(-q/(q+1)^2))*C
      F[a1,1] <- F[a1,1] + A
      F[a2,1] <- F[a2,1] - A
      F[n+1,1] <- F[n+1,1] + A
      F1[a1,1] <- F1[a1,1] + B
      F1[a2,1] <- F1[a2,1] + B
      F1[n+1,1] <- F1[n+1,1] + B
      ##s1 <- s1 - C*df[i,3]*(p-log(q+1))+df[i,4]*(-log(q+1))
    }

    F[,1] <- F[,1]-F[i0,1]
    F1[,1] <- 1/F1[,1]

    ability[,"xn"] <- ability[,"xn"]-F1*F/(1/(j-0.9)+1)

    ##cat(k,ability[21,1],'\n')
    k <- sum(abs(F[,1]))

    j <- j + 1
    if((k<0.00001)||(j>10000)){
      stop <- 1
    }
  }
  ability_CI <- matrix(0,ncol = 2,nrow = (n+1))
  ability_CI[,1] <- ability
  ability_CI[,2] <- F1^0.5
  return(ability_CI)
}

#' Bradely-Terry Model Likelihood function with Lagrangian and Quadratic Penalty
#'
#' @export

B_T_with_Qua <- function(df,ability,i0,theta,n,v,uij){
  n1 <- n*(n-1)
  stop <- 0
  j <- 1
  F <- matrix(0,nrow = (n+1),ncol = 1)
  F1 <- matrix(0,nrow = (n+1),ncol = 1)

  while(stop==0){
    for(i in 1:n){
      F[i,1] <- v*(-sum(ability[-(n+1),1])+n*ability[i,1]+sum(theta[,i])-sum(theta[i,]))+sum(uij[,i])-sum(uij[i,])
      F1[i,1] <- v*(n-1)
    }
    F[(n+1),1] <- 0
    F1[(n+1),1] <- 0
    s <- 0
    for(i in 1:n1){
      a1 <- df[i,1]
      a2 <- df[i,2]
      x <- ability[a1,"xn"]
      y <- ability[a2,"xn"]
      p <- ability["at.home","xn"]+x-y
      q <- exp(p)
      A <- -df[i,3]*(1/(q+1))-df[i,4]*(-q/(q+1))
      B <- -df[i,3]*(-q/(q+1)^2)-df[i,4]*(-q/(q+1)^2)
      F[a1,1] <- F[a1,1] + A
      F[a2,1] <- F[a2,1] - A
      F[n+1,1] <- F[n+1,1] + A
      F1[a1,1] <- F1[a1,1] + B
      F1[a2,1] <- F1[a2,1] + B
      F1[n+1,1] <- F1[n+1,1] + B
      s <- s - df[i,3]*(p-log(q+1))-df[i,4]*(-log(q+1))
    }

    k <- sum(abs(F[,1]))

    F[,1] <- F[,1]-F[i0,1]
    F1[,1] <- 1/F1[,1]

    ##Warm Start
    ability[,"xn"] <- ability[,"xn"]-F1*F/(1/(j-0.9)+1)

    theta1 <- theta
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        theta1[i,j] <- theta[i,j]-ability[i,1]+ability[j,1]
      }
    }
    s <- s + v/2*sum(theta1^2)
    ##cat(s,'\n')

    j <- j + 1
    if((k<0.00001)||(j>10000)){
      stop <- 1
    }
  }
  return(ability)
}

#' Bradely-Terry Model Weighted Likelihood function with Lagrangian and Quadratic Penalty
#'
#' @export

B_T_with_Qua_u <- function(df,ability,i0,theta,n,v,uij,u){
  n1 <- nrow(df)
  stop <- 0
  j <- 1
  F <- matrix(0,nrow = (n+1),ncol = 1)
  F1 <- matrix(0,nrow = (n+1),ncol = 1)

  ability0 <- matrix(0,nrow = n,ncol = 2)

  for(i in 1:n1){
    ability0[df[i,1],1] <- ability0[df[i,1],1] + df[i,3]
    ability0[df[i,2],1] <- ability0[df[i,2],1] + df[i,4]
    ability0[df[i,1],2] <- ability0[df[i,1],2] + 1
    ability0[df[i,2],2] <- ability0[df[i,2],2] + 1
  }

  k1 <- c()
  k0 <- c()
  for(i in 1:n){
    if(ability0[i,1]==ability0[i,2]){
      k1 <- c(k1,i)
    } else if(ability0[i,1]==0){
      k0 <- c(k0,i)
    }
  }

  while(stop==0){
    for(i in 1:n){
      F[i,1] <- v*(-sum(ability[-(n+1),1])+n*ability[i,1]+sum(theta[,i])-sum(theta[i,]))+sum(uij[,i])-sum(uij[i,])
      F1[i,1] <- v*(n-1)
    }
    F[(n+1),1] <- 0
    F1[(n+1),1] <- 0
    for(i in 1:n1){
      a1 <- df[i,1]
      a2 <- df[i,2]
      C <- exp(-u*df[i,5])
      x <- ability[a1,"xn"]
      y <- ability[a2,"xn"]
      p <- ability["at.home","xn"]+x-y
      q <- exp(p)
      A <- (-df[i,3]*(1/(q+1))-df[i,4]*(-q/(q+1)))*C
      B <- (-df[i,3]*(-q/(q+1)^2)-df[i,4]*(-q/(q+1)^2))*C
      F[a1,1] <- F[a1,1] + A
      F[a2,1] <- F[a2,1] - A
      F[n+1,1] <- F[n+1,1] + A
      F1[a1,1] <- F1[a1,1] + B
      F1[a2,1] <- F1[a2,1] + B
      F1[n+1,1] <- F1[n+1,1] + B
    }

    F[,1] <- F[,1]-F[i0,1]
    F1[,1] <- 1/F1[,1]
    if(!identical(k0,NULL)){
      F[k0,1] <- 0
      if(ability0[i0,1]==0){
        ability[k0,1] <- 0
      } else{
        ability[k0,1] <- -1
      }
    }
    if(!identical(k1,NULL)){
      F[k1,1] <- 0
      ability[k1,1] <- 4
    }
    k <- sum(abs(F[,1]))

    ##Warm Start
    ability[,"xn"] <- ability[,"xn"]-F1*F/(1/(j-0.9)+1)

    j <- j + 1
    if((k<0.00001)||(j>10000)){
      stop <- 1
    }
    ##cat(k,'\n')
  }
  return(ability)
}

#' Bradely-Terry Model Weighted Likelihood function with fixed home parameter, Lagrangian and Quadratic Penalty
#'
#' @export

B_T_with_Qua_u_fhome <- function(df,ability,i0,theta,n,v,uij,u,home){
  n1 <- nrow(df)
  n <- nrow(ability)-1
  ability <- as.matrix(ability[1:n,])
  colnames(ability) <- c("xn")

  stop <- 0
  j <- 1
  F <- matrix(0,nrow = n,ncol = 1)
  F1 <- matrix(0,nrow = n,ncol = 1)

  while(stop==0){
    for(i in 1:n){
      F[i,1] <- v*(-sum(ability[-(n+1),1])+n*ability[i,1]+sum(theta[,i])-sum(theta[i,]))+sum(uij[,i])-sum(uij[i,])
      F1[i,1] <- v*(n-1)
    }
    for(i in 1:n1){
      a1 <- df[i,1]
      a2 <- df[i,2]
      C <- exp(-u*df[i,5])
      x <- ability[a1,"xn"]
      y <- ability[a2,"xn"]
      p <- home+x-y
      q <- exp(p)
      A <- (-df[i,3]*(1/(q+1))-df[i,4]*(-q/(q+1)))*C
      B <- (-df[i,3]*(-q/(q+1)^2)-df[i,4]*(-q/(q+1)^2))*C
      F[a1,1] <- F[a1,1] + A
      F[a2,1] <- F[a2,1] - A
      F1[a1,1] <- F1[a1,1] + B
      F1[a2,1] <- F1[a2,1] + B
    }

    k <- sum(abs(F[,1]))

    F[,1] <- F[,1]-F[i0,1]
    F1[,1] <- 1/F1[,1]

    ##Warm Start
    ability[,"xn"] <- ability[,"xn"]-F1*F/(1/(j-0.9)+1)

    j <- j + 1
    if((k<0.00001)||(j>10000)){
      stop <- 1
    }
    ##cat(k,'\n')
  }
  ability <- rbind(ability,matrix(home))
  rownames(ability)[n+1] <- c("at.home")
  return(ability)
}

#' Adapted Lasso's weight determination function
#'
#' @export

B_T_Wij <- function(dataframe){
  ability <- B_T_Quadratic_D(dataframe)
  n <- nrow(ability)-1
  wij <- matrix(0,nrow = n,ncol = n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(abs(ability[i,1]-ability[j,1])<0.01){
        wij[i,j] <- -1
      } else{
        wij[i,j] <- 1/abs(ability[i,1]-ability[j,1])
      }
    }
  }
  k <- max(wij)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(wij[i,j]==-1){
        wij[i,j] <- k
      }
    }
  }
  return(wij)
}

#' Weighted Adapted Lasso's weight determination function
#'
#' @export

B_T_Wij_u <- function(df,ability,i0,u){
  ability <- B_T_Weighted_home(df,ability,u,i0)
  n <- nrow(ability)-1
  wij <- matrix(0,nrow = n,ncol = n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(abs(ability[i,1]-ability[j,1])<0.01){
        wij[i,j] <- -1
      } else{
        wij[i,j] <- 1/abs(ability[i,1]-ability[j,1])
      }
    }
  }
  k <- max(wij)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(wij[i,j]==-1){
        wij[i,j] <- k
      }
    }
  }
  return(wij)
}

#' Weighted Adapted Lasso's weight determination function with fixed home paramter
#'
#' @export

B_T_Wij_u_fhome <- function(df,ability,i0,n,u,home){
  ability <- B_T_Weighted_fhome(df,ability,home,u,i0)
  wij <- matrix(0,nrow = n,ncol = n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(abs(ability[i,1]-ability[j,1])<0.01){
        wij[i,j] <- -1
      } else{
        wij[i,j] <- 1/abs(ability[i,1]-ability[j,1])
      }
    }
  }
  k <- max(wij)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(wij[i,j]==-1){
        wij[i,j] <- k
      }
    }
  }
  return(wij)
}

#' Bradley-Terry Likelihood Function with Lasso Peanlty's optimization given fixed Lagraingian
#' and Quadratic Penalty
#'
#' @export

B_T_Step_1_2 <- function(df,ability,i0,wij,uij,theta,n,v,lambda){
  stop <- 0
  j <- 1
  s1 <- 1000
  while(stop==0){
    ability <- B_T_with_Qua(df,ability,i0,theta,n,v,uij)
    theta <- B_T_Theta(ability,wij,uij,n,v,lambda)
    s <- Likelihood_All(df,ability,theta,n,v,wij,uij,lambda)
    j <- j + 1
    if((abs(s-s1) < 0.0001)||(j>10000)){
      stop <- 1
    } else{
      s1 <- s
    }
    ##cat(s1,'\n')
  }
  return(ability)
}

#' Bradley-Terry Weighted Likelihood Function with Lasso Peanlty's optimization given fixed Lagraingian
#' and Quadratic Penalty
#'
#' @export

B_T_Step_1_2_u <- function(df,ability,i0,wij,uij,theta,n,v,lambda,u){
  stop <- 0
  j <- 1
  s1 <- 1000
  while(stop==0){
    ability <- B_T_with_Qua_u(df,ability,i0,theta,n,v,uij,u)
    theta <- B_T_Theta(ability,wij,uij,n,v,lambda)
    s <- Likelihood_All_u(df,ability,theta,n,v,wij,uij,lambda,u)
    j <- j + 1
    if((abs(s-s1) < 0.0001)||(j>10000)){
      stop <- 1
    } else{
      s1 <- s
    }
    ##cat(s1,'\n')
  }
  return(ability)
}

#' Bradley-Terry Weighted Likelihood Function with fixed home parameter and Lasso Peanlty's optimization given fixed Lagraingian
#' and Quadratic Penalty
#'
#' @export

B_T_Step_1_2_u_fhome <- function(df,ability,i0,wij,uij,theta,n,v,lambda,u,home){
  stop <- 0
  j <- 1
  s1 <- 1000
  while(stop==0){
    ability <- B_T_with_Qua_u_fhome(df,ability,i0,theta,n,v,uij,u,home)
    theta <- B_T_Theta(ability,wij,uij,n,v,lambda)
    s <- Likelihood_All_u(df,ability,theta,n,v,wij,uij,lambda,u)
    j <- j + 1
    if((abs(s-s1) < 0.0001)||(j>10000)){
      stop <- 1
    } else{
      s1 <- s
    }
    ##cat(s1,'\n')
  }
  return(ability)
}

#' Bradley-Terry Likelihood Function with Lasso Peanlty's optimization
#'
#' @param dataframe Dataframe with 4 columns. First column is the home teams
#' Second column is the away teams.
#' Third column is the number of wins of home teams.
#' Fourth column is the number of wins of away teams.
#' @param lambda Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
#' Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.
#' @param wij The weights added to the Lasso Penalty. Can be manually setted or determined using
#' function B_T_Wij
#' @export

B_T_Lasso <- function(dataframe,lambda,wij){
  df1 <- dataframe
  df <- matrix(ncol = ncol(df1),nrow = nrow(df1))
  team <- as.character(unique(df1[,1]))
  n <- length(team)
  for(i in 1:n){
    df[df1==team[i]] <- i
  }
  df[,3] <- df1[,3]
  df[,4] <- df1[,4]
  n1 <- n*(n-1)
  k0 <- 1000
  i0 <- 1
  for(i in 1:n){
    k1 <- sum(df[df[,1]==i,3])-sum(df[df[,1]==i,4])+sum(df[df[,2]==i,4])-sum(df[df[,2]==i,3])
    if(k0>k1){
      i0 <- i
      k0 <- k1
    }
  }

  ##Initialized
  ability <- matrix(0,ncol=1,nrow=(n+1))
  colnames(ability) <- c("xn")
  rownames(ability) <- c(team,"at.home")
  theta <- matrix(0,nrow = n,ncol = n)
  uij <- matrix(0,nrow = n,ncol = n)

  stop <- 0
  j <- 1
  v <- 10
  while(stop==0){
    ability <- B_T_Step_1_2(df,ability,i0,wij,uij,theta,n,v,lambda)
    theta <- B_T_Theta(ability,wij,uij,n,v,lambda)
    uij0 <- B_T_U(uij,ability,theta,v)
    k <- sum(abs(uij0-uij))
    if((k < 0.000001)||(j>100000)){
      stop <- 1
    } else{
      uij <- uij0
      v <- max(uij^2)
    }
    s <- Likelihood_All(df,ability,theta,n,v,wij,uij,lambda)
  }
  return(ability)
}

#' Bradley-Terry Weighted Likelihood Function with Lasso Peanlty's optimization
#'
#' @param dataframe Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @param lambda Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
#' Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.
#' @param wij The weights added to the Lasso Penalty. Can be manually setted or determined using
#' function B_T_Wij
#' @param u The exponential decay rate
#' @return Estimated abilities
#' @export

B_T_Lasso_u <- function(df,ability,i0,lambda,wij,u){
  ##Initialized
  n <- nrow(ability)-1
  theta <- matrix(0,nrow = n,ncol = n)
  uij <- matrix(0,nrow = n,ncol = n)
  ability[,1] <- 0

  stop <- 0
  j <- 1
  v <- 10
  while(stop==0){
    ability <- B_T_Step_1_2_u(df,ability,i0,wij,uij,theta,n,v,lambda,u)
    theta <- B_T_Theta(ability,wij,uij,n,v,lambda)
    uij0 <- B_T_U(uij,ability,theta,v)
    k <- sum(abs(uij0-uij))
    if((k < 0.0001)||(j>100000)){
      stop <- 1
    } else{
      uij <- uij0
      v <- max(uij^2)
    }
    s <- Likelihood_All_u(df,ability,theta,n,v,wij,uij,lambda,u)
  }
  cat(s,'\n')
  return(ability)
}

#' Bradley-Terry Weighted Likelihood Function with fixed home parameter and Lasso Peanlty's optimization
#'
#' @param dataframe Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @param lambda Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
#' Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.
#' @param wij The weights added to the Lasso Penalty. Can be manually setted or determined using
#' function B_T_Wij.
#' @param u The exponential decay rate.
#' @param home Fixed home parameter.
#' @return Estimated abilities
#' @export

B_T_Lasso_u_fhome <- function(df,ability,i0,lambda,wij,u,home){
  ##Initialized
  n <- nrow(ability)-1
  theta <- matrix(0,nrow = n,ncol = n)
  uij <- matrix(0,nrow = n,ncol = n)

  stop <- 0
  j <- 1
  v <- 10
  while(stop==0){
    ability <- B_T_Step_1_2_u_fhome(df,ability,i0,wij,uij,theta,n,v,lambda,u,home)
    theta <- B_T_Theta(ability,wij,uij,n,v,lambda)
    uij0 <- B_T_U(uij,ability,theta,v)
    k <- sum(abs(uij0-uij))
    if((k < 0.000001)||(j>100000)){
      stop <- 1
    } else{
      uij <- uij0
      v <- max(uij^2)
    }
    s <- Likelihood_All_u(df,ability,theta,n,v,wij,uij,lambda,u)
  }
  cat(s,'\n')
  return(ability)
}

#' Bradley-Terry Weighted Likelihood Estimation before a specified time
#' Prediction on one future period is done
#' @export

B_T_Weighted_Lasso_TR <- function(df,ability,i0,date,u,dd){
  n <- nrow(ability)-1
  n1 <- nrow(df)
  n3 <- length(u)

  ##df0 <- df[(date[dd]+1):n1,]
  ##df1 <- df[1:(date[dd]),]
  ##df0[,5] <- df0[,5]-df0[1,5]
  df0 <- df[(date[dd]+1):n1,]
  if(dd==1){
    df1 <- df[1:(date[dd]),]
  } else{
    df1 <- df[(date[dd-1]+1):(date[dd]),]
  }
  df0[,5] <- df0[,5]-df0[1,5]

  ability <- B_T_Weighted_home(df0,ability,0,i0)
  s0 <- Likelihood_Pure(df1,ability,nrow(df1))

  ss0 <- rbind(matrix(ability),matrix(c(s0),ncol=1))

  for(h1 in 1:n3){
    u0 <- u[h1]

    ability <- B_T_Weighted_home(df0,ability,u0,i0)
    s0 <- Likelihood_Pure(df1,ability,nrow(df1))

    ss0 <- cbind(ss0,rbind(matrix(ability),matrix(c(s0),ncol=1)))
  }
  return(ss0)
}

B_T_Weighted_Lasso_TRB <- function(df,ability,i0,date,u,dd){
  n <- nrow(ability)-1
  n1 <- nrow(df)
  n3 <- length(u)

  df0 <- df[(date[dd]+1):n1,]
  if(dd==1){
    df1 <- df[1:(date[dd]),]
  } else{
    df1 <- df[(date[dd-1]+1):(date[dd]),]
  }
  df0[,5] <- df0[,5]-df0[1,5]

  ability <- B_T_Weighted_home(df0,ability,0,i0)
  s0 <- Likelihood_Pure(df1,ability,nrow(df1))

  ss0 <- rbind(matrix(ability),matrix(c(s0),ncol=1))

  for(h1 in 1:n3){
    u0 <- u[h1]

    ability <- B_T_Weighted_home(df0,ability,u0,i0)
    s0 <- Likelihood_Pure(df1,ability,nrow(df1))

    ss0 <- cbind(ss0,rbind(matrix(ability),matrix(c(s0),ncol=1)))
  }
  return(ss0)
}

#' Bradley-Terry Weighted Likelihood function with and Lasso Peanlty's optimization where Wij are setting to be 1
#' Prediction on one future period is done
#'
#' @param dataframe Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @param date A vector which separate all rows of datasets into many parts.
#' @param lambda Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
#' Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.
#' @param u The exponential decay rate.
#' @param dd The index of one value from date. Trained model will be run before this date and prediction will be done
#' on the previous part of dataset.
#' @return The estimated abilities and predicted likelihood
#' @export

B_T_Weighted_Lasso_R <- function(df,ability,i0,date,u,lambda,dd){
  n <- nrow(ability)-1
  n1 <- nrow(df)
  n3 <- length(u)

  df0 <- df[(date[dd]+1):n1,]
  df1 <- df[1:(date[dd]),]
  df0[,5] <- df0[,5]-df0[1,5]

  ##wij <- B_T_Wij_u(df0,ability,i0,n,0)
  wij <- matrix(0,nrow = n,ncol = n)
  for(i in 1:(n-1)){
    wij[1:i,(i+1)] <- 1
  }

  ability <- B_T_Lasso_u(df0,ability,i0,lambda,wij,0)
  k0 <- Lambda_s(ability,wij)
  s0 <- Likelihood_Pure(df1,ability,nrow(df1))
  ss <- rbind(ability,matrix(c(s0,k0),ncol = 1))
  cat("-----------",'\n')

  s1 <- Likelihood_Pure_u(df0,ability,nrow(df0),0)
  for(h1 in 1:n3){
    u0 <- u[h1]

    s2 <- Likelihood_Pure_u(df0,ability,nrow(df0),u0)
    lambda1 <- lambda*s2/s1
    ability <- B_T_Lasso_u(df0,ability,i0,lambda1,wij,u0)
    k1 <- Lambda_s(ability,wij)

    lambda2 <- lambda1*k1/k0
    ability <- B_T_Lasso_u(df0,ability,i0,lambda2,wij,u0)
    k2 <- Lambda_s(ability,wij)

    stop <- 0
    m <- 1
    while(stop==0){
      lambda3 <- (lambda1*k2-lambda2*k1+k0*(lambda2-lambda1))/(k2-k1)
      ability <- B_T_Lasso_u(df0,ability,i0,lambda3,wij,u0)
      k3 <- Lambda_s(ability,wij)
      m <- m+1

      if((abs(k3-k0)<0.01)||(m>15)){
        stop <- 1
      } else{
        if(abs(k0-k1)>abs(k0-k2)){
          k1 <- k2
          lambda1 <- lambda2
        }
        k2 <- k3
        lambda2 <- lambda3
      }
    }
    s0 <- Likelihood_Pure(df1,ability,nrow(df1))
    cat("-----------",'\n')
    ss <- cbind(ss,rbind(ability,matrix(c(s0,k3),ncol=1)))
  }
  return(ss)
}

#' Bradley-Terry Weighted Likelihood function with and Lasso Peanlty's optimization
#' Prediction on one future period is done
#'
#' @param dataframe Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @param date A vector which separate all rows of datasets into many parts.
#' @param lambda Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
#' Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.
#' @param u The exponential decay rate.
#' @param penalty Penalty level = s/max(s)
#' @param dd The index of one value from date. Trained model will be run before this date and prediction will be done
#' on the previous part of dataset.
#' @param uk Exponetial rate between lambda and true penalty level (usually setting to be -0.642)
#' @return The estimated abilities and predicted likelihood
#' @export

B_T_Weighted_Lasso_R_wij <- function(df,ability,i0,date,u,penalty,dd,uk){
  n <- nrow(ability)-1
  n1 <- nrow(df)
  n3 <- length(u)
  k0 <- n*(n-1)/2*penalty
  lambda <- -1/uk*log(k0/(n*(n-1)))
  lambda1 <- lambda

  df0 <- df[(date[dd]+1):n1,]
  if(dd==1){
    df1 <- df[1:(date[dd]),]
  } else{
    df1 <- df[(date[dd-1]+1):(date[dd]),]
  }
  df0[,5] <- df0[,5]-df0[1,5]

  wij <- B_T_Wij_u(df0,ability,i0,0)
  ability <- B_T_Lasso_u(df0,ability,i0,lambda1,wij,0)
  k1 <- Lambda_s(ability,wij)

  lambda2 <- lambda1*k1/k0
  ability <- B_T_Lasso_u(df0,ability,i0,lambda2,wij,0)
  k2 <- Lambda_s(ability,wij)

  stop <- 0
  m <- 1
  while(stop==0){
    if((abs(k1)+abs(k2))<1){
      lambda3 <- 0.5*min(lambda1,lambda2)
    } else{
      lambda3 <- max((lambda1*k2-lambda2*k1+k0*(lambda2-lambda1))/(k2-k1),0)
    }
    ability <- B_T_Lasso_u(df0,ability,i0,lambda3,wij,0)
    k3 <- Lambda_s(ability,wij)
    m <- m+1

    if((abs(k3-k0)<0.01)||(m>15)){
      stop <- 1
    } else{
      if(abs(k0-k1)>abs(k0-k2)){
        k1 <- k2
        lambda1 <- lambda2
      }
      k2 <- k3
      lambda2 <- lambda3
    }
  }

  lambda <- lambda3
  s0 <- Likelihood_Pure(df1,ability,nrow(df1))
  ss <- rbind(ability,matrix(c(s0,k0),ncol = 1))
  cat("-----------",'\n')

  s1 <- Likelihood_Pure_u(df0,ability,nrow(df0),0)
  ability0 <- ability
  for(h1 in 1:n3){
    u0 <- u[h1]

    s2 <- Likelihood_Pure_u(df0,ability0,nrow(df0),u0)
    lambda1 <- lambda*s2/s1

    wij <- B_T_Wij_u(df0,ability,i0,u0)

    ability <- B_T_Lasso_u(df0,ability,i0,lambda1,wij,u0)
    k1 <- Lambda_s(ability,wij)

    lambda2 <- lambda1*k1/k0
    ability <- B_T_Lasso_u(df0,ability,i0,lambda2,wij,u0)
    k2 <- Lambda_s(ability,wij)

    stop <- 0
    m <- 1
    while(stop==0){
      if((abs(k1)+abs(k2))<1){
        lambda3 <- 0.5*min(lambda1,lambda2)
      } else{
        lambda3 <- max((lambda1*k2-lambda2*k1+k0*(lambda2-lambda1))/(k2-k1),0)
      }
      ability <- B_T_Lasso_u(df0,ability,i0,lambda3,wij,u0)
      k3 <- Lambda_s(ability,wij)
      m <- m+1

      if((abs(k3-k0)<0.01)||(m>15)){
        stop <- 1
      } else{
        if(abs(k0-k1)>abs(k0-k2)){
          k1 <- k2
          lambda1 <- lambda2
        }
        k2 <- k3
        lambda2 <- lambda3
      }
    }
    s0 <- Likelihood_Pure(df1,ability,nrow(df1))
    cat("-----------",'\n')
    ss <- cbind(ss,rbind(ability,matrix(c(s0,k3),ncol=1)))
  }
  return(ss)
}

#' (Important) Bradley-Terry Weighted Likelihood Estimation on specifed Lasso Peanlty Level
#'
#' @param dataframe Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @param lambda Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
#' Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.
#' @param u The exponential decay rate.
#' @param penalty Penalty level = s/max(s)
#' @param uk Exponetial rate between lambda and true penalty level (usually setting to be -0.642)
#' @return The estimated abilities
#' @export

B_T_Weighted_Lasso_R_wij_f <- function(df,ability,i0,u,penalty,uk){
  n <- nrow(ability)-1
  n1 <- nrow(df)
  n3 <- length(u)
  k0 <- n*(n-1)/2*penalty
  lambda <- -1/uk*log(k0/(n*(n-1)))
  lambda1 <- lambda

  df0 <- df
  df0[,5] <- df0[,5]-df0[1,5]

  wij <- B_T_Wij_u(df0,ability,i0,0)
  ability <- B_T_Lasso_u(df0,ability,i0,lambda1,wij,0)
  k1 <- Lambda_s(ability,wij)

  lambda2 <- lambda1*k1/k0
  ability <- B_T_Lasso_u(df0,ability,i0,lambda2,wij,0)
  k2 <- Lambda_s(ability,wij)

  stop <- 0
  m <- 1
  while(stop==0){
    if((abs(k1)+abs(k2))<1){
      lambda3 <- 0.5*min(lambda1,lambda2)
    } else{
      lambda3 <- max((lambda1*k2-lambda2*k1+k0*(lambda2-lambda1))/(k2-k1),0)
    }
    ability <- B_T_Lasso_u(df0,ability,i0,lambda3,wij,0)
    k3 <- Lambda_s(ability,wij)
    m <- m+1

    if(lambda3==0&&k3<k0){
      stop <- 1
      lambda3 <- 0.01
    }

    if((abs(k3-k0)<0.01)||(m>15)){
      stop <- 1
    } else{
      if(abs(k0-k1)>abs(k0-k2)){
        k1 <- k2
        lambda1 <- lambda2
      }
      k2 <- k3
      lambda2 <- lambda3
    }
  }

  lambda <- lambda3
  ss <- rbind(ability,matrix(c(k3),ncol = 1))
  cat("-----------",'\n')

  s1 <- Likelihood_Pure_u(df0,ability,nrow(df0),0)
  ability0 <- ability
  for(h1 in 1:n3){
    u0 <- u[h1]

    s2 <- Likelihood_Pure_u(df0,ability0,nrow(df0),u0)
    lambda1 <- lambda*s2/s1

    wij <- B_T_Wij_u(df0,ability,i0,u0)

    ability <- B_T_Lasso_u(df0,ability,i0,lambda1,wij,u0)
    k1 <- Lambda_s(ability,wij)

    lambda2 <- lambda1*k1/k0
    ability <- B_T_Lasso_u(df0,ability,i0,lambda2,wij,u0)
    k2 <- Lambda_s(ability,wij)

    stop <- 0
    m <- 1
    while(stop==0){
      if((abs(k1)+abs(k2))<1){
        lambda3 <- 0.5*min(lambda1,lambda2)
      } else{
        lambda3 <- max((lambda1*k2-lambda2*k1+k0*(lambda2-lambda1))/(k2-k1),0)
      }
      ability <- B_T_Lasso_u(df0,ability,i0,lambda3,wij,u0)
      k3 <- Lambda_s(ability,wij)
      m <- m+1

      if(lambda3==0&&k3<k0){
        stop <- 1
        lambda3 <- 0.01
      }

      if((abs(k3-k0)<0.01)||(m>15)){
        stop <- 1
      } else{
        if(abs(k0-k1)>abs(k0-k2)){
          k1 <- k2
          lambda1 <- lambda2
        }
        k2 <- k3
        lambda2 <- lambda3
      }
    }
    cat("-----------",'\n')
    ss <- cbind(ss,rbind(ability,matrix(c(k3),ncol=1)))
  }
  return(ss)
}

#' Bradley-Terry Weighted Likelihood Estimation with fixed home parameter and Lasso Peanlty's optimization
#'
#' @param dataframe Dataframe with 5 columns. First column is the index of the home teams
#' (use numbers to denote teams).
#' Second column is the index of the away teams.
#' Third column is the number of wins of home teams (usually to be 0/1).
#' Fourth column is the number of wins of away teams (usually to be 0/1).
#' Fifth column is the scalar of time when the match is played until now (Time lag).
#' @param ability A column vector of teams ability, the last row is the home parameter whose rowname must be "at.home".
#' @param i0 A teams index whose ability will be fixed as 0 (usually the team loss most).
#' @param lambda Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
#' Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.
#' @param u The exponential decay rate.
#' @param penalty Penalty level = s/max(s)
#' @param uk Exponetial rate between lambda and true penalty level (usually setting to be -0.642)
#' @param home Fixed home ability
#' @return The estimated abilities
#' @export

B_T_Weighted_Lasso_R_wij_fhome <- function(df,ability,i0,u,penalty,uk,home){
  n <- nrow(ability)-1
  n1 <- nrow(df)
  n3 <- length(u)
  k0 <- n*(n-1)/2*penalty
  lambda <- -1/uk*log(k0/(n*(n-1)))
  lambda1 <- lambda

  df0 <- df
  df0[,5] <- df0[,5]-df0[1,5]

  wij <- B_T_Wij_u(df0,ability,i0,0)
  ability <- B_T_Lasso_u_fhome(df0,ability,i0,lambda1,wij,0,home)
  k1 <- Lambda_s(ability,wij)

  lambda2 <- lambda1*k1/k0
  ability <- B_T_Lasso_u_fhome(df0,ability,i0,lambda2,wij,0,home)
  k2 <- Lambda_s(ability,wij)

  stop <- 0
  m <- 1
  lambda3 <- 0
  while(stop==0){
    lambda3_1 <- lambda3
    if((abs(k1)+abs(k2))<0.1){
      lambda3 <- 0.5*min(lambda1,lambda2)
    } else{
      lambda3 <- max((lambda1*k2-lambda2*k1+k0*(lambda2-lambda1))/(k2-k1),0)
    }
    ability <- B_T_Lasso_u_fhome(df0,ability,i0,lambda3,wij,0,home)
    k3 <- Lambda_s(ability,wij)
    m <- m+1

    if((abs(k3-k0)<0.01)||(m>15)||(abs(lambda3_1-lambda3)<1e-20)){
      stop <- 1
    } else{
      if(abs(k0-k1)>abs(k0-k2)){
        k1 <- k2
        lambda1 <- lambda2
      }
      k2 <- k3
      lambda2 <- lambda3
    }
  }

  lambda <- lambda3
  ss <- rbind(ability,matrix(c(k3),ncol = 1))
  cat("-----------",'\n')

  s1 <- Likelihood_Pure_u(df0,ability,nrow(df0),0)
  ability0 <- ability
  for(h1 in 1:n3){
    u0 <- u[h1]

    s2 <- Likelihood_Pure_u(df0,ability0,nrow(df0),u0)
    lambda1 <- lambda*s2/s1

    wij <- B_T_Wij_u(df0,ability,i0,u0)

    ability <- B_T_Lasso_u_fhome(df0,ability,i0,lambda1,wij,u0,home)
    k1 <- Lambda_s(ability,wij)

    lambda2 <- lambda1*k1/k0
    ability <- B_T_Lasso_u_fhome(df0,ability,i0,lambda2,wij,u0,home)
    k2 <- Lambda_s(ability,wij)

    stop <- 0
    m <- 1
    lambda3 <- 0
    while(stop==0){
      lambda3_1 <- lambda3
      if((abs(k1)+abs(k2))<0.1){
        lambda3 <- 0.5*min(lambda1,lambda2)
      } else{
        lambda3 <- max((lambda1*k2-lambda2*k1+k0*(lambda2-lambda1))/(k2-k1),0)
      }
      ability <- B_T_Lasso_u_fhome(df0,ability,i0,lambda3,wij,u0,home)
      k3 <- Lambda_s(ability,wij)
      m <- m+1

      if((abs(k3-k0)<0.01)||(m>15)||(abs(lambda3_1-lambda3)<1e-20)){
        stop <- 1
      } else{
        if(abs(k0-k1)>abs(k0-k2)){
          k1 <- k2
          lambda1 <- lambda2
        }
        k2 <- k3
        lambda2 <- lambda3
      }
    }
    cat("-----------",'\n')
    ss <- cbind(ss,rbind(ability,matrix(c(k3),ncol=1)))
  }
  return(ss)
}

#' Theta update Function
#' @export

B_T_Theta <- function(ability,wij,uij,n,v,lambda){
  theta <- matrix(0,nrow = n,ncol = n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      theta0 <- ability[i,1]-ability[j,1]-uij[i,j]/v
      theta[i,j] <- sign(theta0)*max(abs(theta0)-lambda*wij[i,j]/v,0)
    }
  }
  return(theta)
}

#' Likelihood computing Function on Bradley-Terry Likelihood only
#' @export

Likelihood_Pure <- function(df,ability,n){
  n1 <- nrow(df)
  s <- 0
  for(i in 1:n1){
    a1 <- df[i,1]
    a2 <- df[i,2]
    x <- ability[a1,"xn"]
    y <- ability[a2,"xn"]
    p <- ability["at.home","xn"]+x-y
    q <- exp(p)
    s <- s - df[i,3]*(p-log(q+1))-df[i,4]*(-log(q+1))
  }
  return(s)
}

#' Likelihood computing Function on weighted Bradley-Terry Likelihood only
#' @export

Likelihood_Pure_u <- function(df,ability,n,u){
  n1 <- nrow(df)
  s <- 0
  for(i in 1:n1){
    a1 <- df[i,1]
    a2 <- df[i,2]
    C <- exp(-u*df[i,5])
    x <- ability[a1,"xn"]
    y <- ability[a2,"xn"]
    p <- ability["at.home","xn"]+x-y
    q <- exp(p)
    s <- s - (df[i,3]*(p-log(q+1))+df[i,4]*(-log(q+1)))*C
  }
  return(s)
}

#' Likelihood computing Function on Bradley-Terry Likelihood and Lasso penalty Part
#' @export

Likelihood <- function(df,ability,wij,lambda,theta,n){
  n1 <- n*(n-1)
  s <- 0
  for(i in 1:n1){
    a1 <- df[i,1]
    a2 <- df[i,2]
    x <- ability[a1,"xn"]
    y <- ability[a2,"xn"]
    p <- ability["at.home","xn"]+x-y
    q <- exp(p)
    s <- s - df[i,3]*(p-log(q+1))-df[i,4]*(-log(q+1))
  }
  s <- s + lambda*sum(abs(theta)*wij)
  return(s)
}

#' Likelihood computing Function on Bradley-Terry Likelihood and Quadratic penalty Part
#' @export

Likelihood_Qua <- function(df,ability,theta,n,v){
  n1 <- n*(n-1)
  s <- 0
  for(i in 1:n1){
    a1 <- df[i,1]
    a2 <- df[i,2]
    x <- ability[a1,"xn"]
    y <- ability[a2,"xn"]
    p <- ability["at.home","xn"]+x-y
    q <- exp(p)
    s <- s - df[i,3]*(p-log(q+1))-df[i,4]*(-log(q+1))
  }
  theta1 <- theta
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      theta1[i,j] <- theta[i,j]-ability[i,1]+ability[j,1]
    }
  }
  s <- s + v/2*sum(theta1^2)
  return(s)
}

#' Likelihood computing Function on Bradley-Terry Likelihood and Lasso and Quadratic penalty Part
#' @export

Likelihood_Qua_lambda <- function(df,ability,theta,n,v,wij,lambda){
  n1 <- n*(n-1)
  s <- 0
  for(i in 1:n1){
    a1 <- df[i,1]
    a2 <- df[i,2]
    x <- ability[a1,"xn"]
    y <- ability[a2,"xn"]
    p <- ability["at.home","xn"]+x-y
    q <- exp(p)
    s <- s - df[i,3]*(p-log(q+1))-df[i,4]*(-log(q+1))
  }
  theta1 <- theta
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      theta1[i,j] <- theta[i,j]-ability[i,1]+ability[j,1]
    }
  }
  s <- s + v/2*sum(theta1^2)
  s <- s + lambda*sum(abs(theta)*wij)
  return(s)
}

#' Likelihood computing Function on Bradley-Terry Likelihood and all penalty Part
#' @export

Likelihood_All <- function(df,ability,theta,n,v,wij,uij,lambda){
  n1 <- nrow(df)
  s <- 0
  for(i in 1:n1){
    a1 <- df[i,1]
    a2 <- df[i,2]
    x <- ability[a1,"xn"]
    y <- ability[a2,"xn"]
    p <- ability["at.home","xn"]+x-y
    q <- exp(p)
    s <- s - df[i,3]*(p-log(q+1))-df[i,4]*(-log(q+1))
  }
  theta1 <- theta
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      theta1[i,j] <- theta[i,j]-ability[i,1]+ability[j,1]
    }
  }
  s <- s + v/2*sum(theta1^2) + sum(uij*theta1)
  s <- s + lambda*sum(abs(theta)*wij)
  return(s)
}

#' Likelihood computing Function on Weighted Bradley-Terry Likelihood and all penalty Part
#' @export

Likelihood_All_u <- function(df,ability,theta,n,v,wij,uij,lambda,u){
  n1 <- nrow(df)
  s <- 0
  for(i in 1:n1){
    a1 <- df[i,1]
    a2 <- df[i,2]
    C <- exp(-u*df[i,5])
    x <- ability[a1,"xn"]
    y <- ability[a2,"xn"]
    p <- ability["at.home","xn"]+x-y
    q <- exp(p)
    s <- s - (df[i,3]*(p-log(q+1))+df[i,4]*(-log(q+1)))*C
  }
  theta1 <- theta
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      theta1[i,j] <- theta[i,j]-ability[i,1]+ability[j,1]
    }
  }
  s <- s + v/2*sum(theta1^2) + sum(uij*theta1)
  s <- s + lambda*sum(abs(theta)*wij)
  return(s)
}

#' Penalty Determine Function
#' @export

Lambda_s <- function(ability,wij){
  n <- nrow(ability)-1
  s <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      s <- s + abs(ability[i,1]-ability[j,1])*wij[i,j]
    }
  }
  return(s)
}

#' Lagrangian Update Function
#' @export

B_T_U <- function(uij,ability,theta,v){
  n <- nrow(ability)-1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      uij[i,j] <- uij[i,j]+v*(theta[i,j]-ability[i,1]+ability[j,1])
    }
  }
  return(uij)
}

#' MSE of linearly changed ability with two teams involved and no home effect
#'
#' @param a Exponetial decay rate.
#' @param n Number of matched.
#' @param a1 The ability score at the beginning.
#' @param k1 Slope of linearly changed ability
#' @param p Lasso penalty level = s/max(s)
#' @return Bias, Variance and MSE
#' @export

MSE_v2 <- function(a,n,a1,k1,p){
  n1 <- 2^n
  s0 <- 0
  s1 <- 0
  pk <- c()
  ak <- c()
  for(i in 1:n){
    pk <- c(pk,(a1+(1-(i-1)/n)*k1)/(a1+(1-(i-1)/n)*k1+1))
    ak <- c(ak,exp(-(i-1)/n*a))
  }
  B <- sum(ak)
  for(i in 0:(n1-1)){
    k <- as.integer(intToBits(i))[1:n]
    A <- sum(k * ak)
    s0 <- s0 + ((2*p-1)*A+(1-p)*B)/((2*p-2)*A+(2-p)*B)*prod(pk*k+(1-pk)*(1-k))
    s1 <- s1 + (((2*p-1)*A+(1-p)*B)/((2*p-2)*A+(2-p)*B))^2*prod(pk*k+(1-pk)*(1-k))
  }
  var <- s1 - s0^2
  bia <- (s0 - (a1+k1)/(a1+k1+1))^2
  mse <- bia + var
  M <- matrix(c(bia, var, mse))
  rownames(M) <- c("bias","var","mse")
  return(M)
}

#' MSE of changed ability follows step function with two teams involved and no home effect
#'
#' @param a Exponetial decay rate.
#' @param n Number of matched.
#' @param a1 The ability score at the beginning.
#' @param k1 Gap of step changed ability
#' @param p Lasso penalty level = s/max(s)
#' @param m Proportion of most recent unchanged ability
#' @return Bias, Variance and MSE
#' @export

MSE_v4 <- function(a,n,a1,k1,p,m){
  n1 <- 2^n
  s0 <- 0
  s1 <- 0
  pk <- c()
  ak <- c()
  for(i in 1:n){
    ak <- c(ak,exp(-(i-1)/n*a))
  }
  for(i in 1:m){
    pk <- c(pk,(a1+k1)/(a1+k1+1))
  }
  for(i in (m+1):n){
    pk <- c(pk,(a1)/(a1+1))
  }
  B <- sum(ak)
  for(i in 0:(n1-1)){
    k <- as.integer(intToBits(i))[1:n]
    A <- sum(k * ak)
    s0 <- s0 + ((2*p-1)*A+(1-p)*B)/((2*p-2)*A+(2-p)*B)*prod(pk*k+(1-pk)*(1-k))
    s1 <- s1 + (((2*p-1)*A+(1-p)*B)/((2*p-2)*A+(2-p)*B))^2*prod(pk*k+(1-pk)*(1-k))
  }
  var <- s1 - s0^2
  bia <- (s0 - (a1+k1)/(a1+k1+1))^2
  mse <- bia + var
  M <- matrix(c(bia, var, mse))
  rownames(M) <- c("bias","var","mse")
  return(M)
}

#' Time Table construction for two teams
#'
#' @export

Time_table0 <- function(x1,x2){
  n <- length(x1)
  d <- seq(0,1-1/n,1/n)
  w <- rep(0,n)
  l <- rep(0,n)
  M <- matrix(c(x1,x2,w,l,d),ncol = 5)
  for(i in 1:n){
    M[i,3:4] <- Pwl(0,M[i,1],M[i,2])
  }
  M[,1] <- 1
  M[,2] <- 2
  return(M)
}

#' Time Table construction for two teams with home effect
#'
#' @export

Time_table2 <- function(n,x1,x2,home){
  K <- matrix(0,ncol = 5,nrow = 2*n)
  for(i in 1:(2*n)){
    if(i/2==floor(i/2)){
      K[i,1] <- 2
      K[i,2] <- 1
      K[i,3:4] <- Pwl(home, x2, x1)
    } else{
      K[i,1] <- 1
      K[i,2] <- 2
      K[i,3:4] <- Pwl(home, x1, x2)
    }
    K[i,5] <- (i-1)/(2*n)
  }
  return(K)
}

#' Time Table construction for four teams with home effect
#'
#' @export

Time_table4 <- function(n,x1,x2,x3,x4,home){
  T <- c(x1,x2,x3,x4)
  K <- matrix(0,ncol = 5,nrow = 12)
  K[1,1:2] <- c(1,3)
  K[2,1:2] <- c(2,4)
  K[3,1:2] <- c(1,4)
  K[4,1:2] <- c(3,2)
  K[5,1:2] <- c(1,2)
  K[6,1:2] <- c(4,3)
  K[7:12,1] <- K[1:6,2]
  K[7:12,2] <- K[1:6,1]
  K1 <- K
  for(i in 1:n){
    K1 <- rbind(K1,K)
  }
  n1 <- nrow(K1)
  for(i in 1:n1){
    K1[i,3:4] <- Pwl(home,T[K1[i,1]],T[K1[i,2]])
    K1[i,5] <- floor((i-1)/2)*2/n1
  }
  return(K1)
}

#' Random Bernulli Trials Generating Functions with home effect
#'
#' @export

Pwl <- function(home,a,b){
  p <- exp(home+a-b)/(1+exp(home+a-b))
  pa <- runif(1,0,1)
  if(pa < p){
    return(matrix(c(1,0),nrow = 1))
  }else{
    return(matrix(c(0,1),nrow = 1))
  }
}

#' Multi-season Data Merge
#'
#' @export

PL_merge <- function(df1,df2,dt){
  for(i in 1:3){
    x1 <- which(df2[,1]==dt[i,2])
    df2[x1,1] <- dt[i,1]
    x2 <- which(df2[,2]==dt[i,2])
    df2[x2,2] <- dt[i,1]
  }
  team <- as.character(unique(df1[,1]))
  df <- matrix(ncol = ncol(df1),nrow = (2*nrow(df1)))
  n <- length(team)

  df0 <- rbind(df1,df2)
  for(i in 1:n){
    df[df0[,1:2]==team[i]] <- i
  }
  df[,3] <- df0[,3]
  df[,4] <- df0[,4]
  df[,5] <- df0[,5]

  x1 <- df[1,5]
  df[,5] <- x1-df[,5]
  df <- df[order(df[,5]),]
  return(df)
}

#' Round clarification Function
#'
#' @export

date_c <- function(df){
  date <- c(0)
  date_t <- as.data.frame(table(df[,5]))
  i <- 1
  k <- 0
  while(i<nrow(date_t)){
    a <- date_t[i,2]
    b <- date_t[(i+1),2]
    if((a+b)<(3/4*n)){
      if((i+2)<nrow(date_t)+1){
        c <- date_t[(i+2),2]
        if((a+b+c)<(3/4*n)){
          k <- k+a+b+c
          date <- c(date,k)
          i <- i + 3
        }else{
          k <- k+a+b
          date <- c(date,k)
          i <- i + 2
        }
      } else{
        i <- i + 2
      }
    }else{
      k <- k+a
      date <- c(date,k)
      i <- i + 1
    }
  }
  date <- date-1
  date[1] <- 0
  date <- c(date,nrow(df))
  date <- date[-1]
  return(date)
}

#' Environment Loading Function
#'
#' @export

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}
