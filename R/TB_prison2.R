TB_prison <- function (t, state, parameters) {
  with(as.list(c(state, parameters)),             
       {
         N       <- U + Lf+ Ls + I + R                 # Total population
         Np      <- Up + Lfp + Lsp + Ip + Rp            # Total population
         Nx      <- Ux + Lfx+ Lsp + Ix + Rx            # Total population
         births <- (I+Ip+Ix) * mu_tb + (N+Np+Nx) * mu  # Births (for stable population)
         lambda   <- beta * (I+Ix)/(N+Nx)        # Force of Infection coomunirt
         lambda_p <- beta_p * Ip/Np              # Force of Infection prisons
         
         
         
         ##### Community transmission
         # Uninfected 
         dU <- births - U * (lambda + mu + r_incar)                              
         
         # Latent fast 
         dLf <- U * lambda + R * lambda * imm - Lf * (mu + break_fast + r_slow + r_incar) 
         
         # Latent slow
         dLs <- Lf * r_slow  - Ls * (mu + break_slow + r_incar) 
         
         # Active TB
         dI <- Lf * break_fast + Ls * break_slow + R * relapse -
           I * (mu + mu_tb + selfcure + r_incar)
         
         # Recovered
         dR <- I * selfcure - R * (relapse + mu + r_incar)   
         
         
         ###### Ex-prison transmission
         # Uninfected 
         dUx <- Up * r_release - Ux * (lambda + mu + r_reincar)                              
         
         # Latent fast 
         dLfx <- Lfp  * r_release + Ux * lambda + Rx * lambda * imm - 
           Lfx * (mu + break_fast_x + r_slow + r_reincar) 
         
         # Latent slow
         dLsx <- Lsp  * r_release + Lfx * r_slow  - Lsx * (mu + break_slow_x + r_reincar) 
         
         # Active TB
         dIx <- Ip * r_release + Lfx * break_fast_x + Lsx * break_slow_x + Rx * relapse -
           Ix * (mu + mu_tb + selfcure_x + r_reincar)
         
         # Recovered
         dRx <- Rp * r_release + Ix * selfcure_x - Rx * (relapse + mu + r_reincar)     
         
         
         ## Prison Transmission
         
         dUp <- U * r_incar + Up * r_reincar - Up * (lambda_p + mu + r_release)                              
         
         # Latent fast 
         dLfp <- Lf * r_incar + Lfx  * r_reincar + Up * lambda_p + Rp * lambda_p * imm - 
           Lfp * (mu + break_fast_p + r_slow + r_release) 
         
         # Latent slow
         dLsp <- Ls * r_incar + Lsx * r_reincar + Lfp * r_slow  - 
           Lsp * (mu + break_slow_p + r_release) 
         
         # Active TB
         dIp <- I * r_incar + Ix * r_reincar + Lfp * break_fast_p + 
           Lsp * break_slow_p + Rp * relapse - Ip * (mu + mu_tb + selfcure_p + r_release)
         
         # Recovered
         dRp <- R * r_incar + Rx * r_reincar + Ip * selfcure_p - Rp * (relapse + mu + r_release)  
         
         
         # Model outcomes
         
         # Community 
         dIncidence <- Lf * break_fast + Ls * break_slow + R * relapse
         
         dIrecent   <- Lf * break_fast
         
         dIremote   <- Ls * break_slow + R * relapse
         
         # Ex Prison 
         dIncidence_x <- Ip * r_release + Lfx * break_fast_x + Lsx * break_slow_x + Rx * relapse
         
         dIrecent_x   <- Lfx * break_fast_x
         
         dIremote_x   <- Lsx * break_slow_x + Rx * relapse
         
         
         # Prison 
         dIncidence_p <-   I * r_incar + Ix * r_reincar + Lfp * break_fast_p + 
           Lsp * break_slow_p + Rp * relapse 
         
         dIrecent_p   <- Lfp * break_fast_p
         
         dIremote_p   <- Lsp * break_slow_p + Rp * relapse 
         
         
         
         # ::::::::::::::::::::::::::::::::::::::::::::
         
         # wrap-up 
         dx <- c(dU, dLf, dLs, dI, dR,
                 dUx, dLfx, dLsx, dIx, dRx,
                 dUp, dLfp,dLsp, dIp, dRp,
                 dIncidence, dIrecent, dIremote,
                 dIncidence_x, dIrecent_x, dIremote_x,
                 dIncidence_p, dIrecent_p, dIremote_p)
         list(dx)
       }
  )
}
