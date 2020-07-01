      module integrator_module

c***********************************************************************
c     
c     dl_poly module for selecting verlet integration schemes
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      use error_module
      use lf_motion_module
      use lf_rotation1_module
      use lf_rotation2_module
      use pmf_module
      use temp_scalers_module
      use vv_motion_module
      use vv_rotation1_module
      use vv_rotation2_module
      
      contains
      
      subroutine lf_integrate
     x  (lcnb,lshmov,lnfic,idnode,mxnode,imcon,natms,nstep,ngrp,
     x  keyens,nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nofic,
     x  tstep,engke,engrot,tolnce,quattol,vircon,vircom,virtot,
     x  temp,press,volm,sigma,taut,taup,chit,chip,consv,conint,
     x  elrc,virlrc,virpmf)

c***********************************************************************
c     
c     dl_poly subroutine for selecting the integration algorithm
c     to solve the the equations of motion. based on the leapfrog
c     verlet algorithm
c     
c     copyright - daresbury laboratory
c     author    - w. smith december 2005
c     
c***********************************************************************

      implicit none

      logical safe,safep,safeq,lcnb,lshmov,lnfic
      integer idnode,mxnode,imcon,natms,ngrp,keyens,nscons,nofic
      integer ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nstep
      real(8) tstep,engke,engrot,tolnce,quattol,vircon,vircom
      real(8) virtot,temp,press,volm,sigma,taut,taup,chit,chip
      real(8) consv,conint,elrc,virlrc,virpmf
      
      safe=.true.
      safeq=.true.
      safep=.true.
      
      if(ngrp.eq.0) then
        
        if(keyens.eq.0) then

c     verlet leapfrog 

          call nve_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,tolnce,tstep,vircon)
          
        else if(keyens.eq.1) then

c     Evans Gaussian Temperature constraints
          
          call nvt_e1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,tolnce,tstep,vircon)
          
        else if(keyens.eq.2) then

c     Berendsen thermostat
          
          call nvt_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,taut,sigma,tolnce,tstep,vircon)
          
        else if(keyens.eq.3) then

c     Nose-Hoover thermostat
          
          call nvt_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      chit,consv,conint,engke,taut,sigma,tolnce,tstep,vircon)
          
        elseif(keyens.eq.4) then

c     Berendsen thermostat and isotropic barostat

          call npt_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,elrc,engke,virlrc,press,taup,taut,sigma,tolnce,
     x      tstep,virtot,vircon,volm)

        else if(keyens.eq.5) then

c     Nose-Hoover thermostat and isotropic barostat 

          call npt_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,chip,chit,conint,consv,elrc,engke,virlrc,press,
     x      taup,taut,sigma,temp,tolnce,tstep,virtot,vircon,volm)

        else if(keyens.eq.6) then

c     Berendsen thermostat and barostat (cell shape varying)

          call nst_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,mode,elrc,engke,virlrc,press,taup,taut,sigma,
     x      tolnce,tstep,vircon,volm)

        else if(keyens.eq.7) then

c     Nose-Hoover thermostat and barostat (cell shape varying)
          
          call nst_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,mode,chit,conint,consv,elrc,engke,virlrc,press,
     x      taup,taut,sigma,temp,tolnce,tstep,vircon,volm)

        elseif(keyens.eq.8) then

c     Potential of mean force in NVE

            call pmflf
     x        (safe,safep,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        ntcons,nspmf,ntpmf,engke,tolnce,tstep,vircon,virpmf)

        endif

      elseif(ngrp.gt.0) then

c     apply rigid body equations of motion
        
        if(keyens.eq.0) then
          
          if(.not.lcnb) then

            call nveq_1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x        vircon)

          else

            call nveq_2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x        vircon)

          endif

        elseif(keyens.eq.1) then

c     invalid option

          call error(idnode,430)
          
        elseif(keyens.eq.2) then
          
          if(.not.lcnb) then

            call nvtq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x        tstep,vircom,vircon)

          else

            call nvtq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x        tstep,vircom,vircon)
          
          endif

        elseif(keyens.eq.3) then
          
          if(.not.lcnb) then 

            call nvtq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,chit,consv,conint,engke,engrot,quattol,
     x        sigma,taut,tolnce,tstep,vircom,vircon)

          else

            call nvtq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,conint,consv,chit,engke,engrot,quattol,
     x        sigma,taut,tolnce,tstep,vircom,vircon)

          endif
            
        elseif(keyens.eq.4) then

          if(.not.lcnb) then

            call nptq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,virtot,vircom,
     x        vircon,volm)
          
          else

            call nptq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,
     x        virtot,volm)

          endif

        elseif(keyens.eq.5) then
          
          if(.not.lcnb) then 

            call nptq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,virtot,vircom,vircon,volm)

          else

            call nptq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,vircom,vircon,virtot,volm)

          endif
            
        elseif(keyens.eq.6) then

          if(.not.lcnb) then

            call nstq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm)

          else

            call nstq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm)

          endif

        elseif(keyens.eq.7) then

          if(.not.lcnb) then

            call nstq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,chit,conint,consv,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,vircom,vircon,volm)

          else

            call nstq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,chit,conint,consv,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,vircom,vircon,volm)
            
          endif

        else

c     invalid option

          call error(idnode,430)

        endif

      endif

c    check on convergence of pmf-shake

      if(ntpmf.gt.0) then

        if(mxnode.gt.1) call gstate(safep)
        if(.not.safep) call error(idnode,438)

      endif    

c    check on convergence of shake

      if(ntcons.gt.0) then

        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe) call error(idnode,105)

      endif    

c     check on convergence of quaternion algorithm

      if(ngrp.gt.0) then

        if(mxnode.gt.1) call gstate(safeq)
        if(.not.safeq) call error(idnode,321)

      endif

c     eliminate "flying ice cube" in long simulations (Berendsen)
      
      if(lnfic.and.(keyens.eq.2.or.keyens.eq.4.or.keyens.eq.6))then
        
        if(mod(nstep,nofic).eq.0)then
          
          call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
          
        endif
        
      endif
      
      return
      end subroutine lf_integrate

      subroutine vv_integrate
     x  (lcnb,lshmov,lnfic,isw,idnode,mxnode,imcon,natms,nstep,ngrp,
     x  keyens,nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nofic,
     x  ntshl,keyshl,tstep,engke,engrot,tolnce,vircon,vircom,virtot,
     x  temp,press,volm,sigma,taut,taup,chit,chip,consv,conint,elrc,
     x  virlrc,virpmf,chit_shl,sigma_shl)

c***********************************************************************
c     
c     dl_poly subroutine for selecting the integration algorithm
c     to solve the the equations of motion. based on the velocity
c     verlet algorithm
c     
c     copyright - daresbury laboratory
c     author    - w. smith february 2005
c     
c***********************************************************************

      implicit none

      logical safe,safep,lcnb,lshmov,lnfic
      integer isw,idnode,mxnode,imcon,natms,ngrp,keyens,nscons
      integer ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nstep,nofic
      integer ntshl,keyshl
      real(8) tstep,engke,engrot,tolnce,vircon,vircom
      real(8) virtot,temp,press,volm,sigma,taut,taup,chit,chip
      real(8) consv,conint,elrc,virlrc,virpmf,chit_shl,sigma_shl
      
      if(ngrp.eq.0) then
        
        if(keyens.eq.0) then

c     verlet leapfrog 

          call nvevv_1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,tstep,engke,tolnce,vircon)
          
        else if(keyens.eq.1) then

c     Evans Gaussian Temperature constraints
          
          call nvtvv_e1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,tstep,engke,tolnce,vircon)
          
        else if(keyens.eq.2) then

c     Berendsen thermostat
          
          call nvtvv_b1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,tstep,taut,sigma,engke,tolnce,vircon)
          
        else if(keyens.eq.3) then

c     Nose-Hoover thermostat
          
          call nvtvv_h1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntshl,keyshl,tstep,taut,sigma,chit,consv,
     x      conint,engke,tolnce,vircon,chit_shl,sigma_shl)
          
        elseif(keyens.eq.4) then

c     Berendsen thermostat and isotropic barostat

          call nptvv_b1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,tstep,taut,taup,sigma,engke,press,elrc,
     x      virlrc,tolnce,virtot,vircon,volm)

        else if(keyens.eq.5) then

c     Nose-Hoover thermostat and isotropic barostat 

          call nptvv_h1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,temp,
     x      chip,chit,consv,conint,engke,elrc,tolnce,vircon,
     x      virtot,virlrc,volm,press,chit_shl,sigma_shl)

        else if(keyens.eq.6) then

c     Berendsen thermostat and barostat (cell shape varying)

          call nstvv_b1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,mode,tstep,taut,taup,sigma,engke,press,
     x      elrc,virlrc,tolnce,vircon,volm)

        else if(keyens.eq.7) then

c     Nose-Hoover thermostat and barostat (cell shape varying)
          
          call nstvv_h1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,mode,ntshl,keyshl,tstep,taut,taup,sigma,
     x      temp,chit,consv,conint,engke,elrc,tolnce,vircon,
     x      virlrc,volm,press,chit_shl,sigma_shl)

        elseif(keyens.eq.8) then

c     Potential of mean force in NVE

          call pmfvv
     x      (safe,safep,lshmov,isw,idnode,mxnode,imcon,natms,nscons,
     x      ntcons,nspmf,ntpmf,engke,tolnce,tstep,vircon,virpmf)

        endif

      elseif(ngrp.gt.0) then

c     apply rigid body equations of motion
        
        if(keyens.eq.0) then
          
          if(.not.lcnb) then

            call nveqvv_1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,tolnce,tstep,vircom,vircon)

          else

            call nveqvv_2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,tolnce,tstep,vircom,vircon)
            
          endif

        elseif(keyens.eq.1) then

c     invalid option

          call error(idnode,430)
          
        elseif(keyens.eq.2) then
          
          if(.not.lcnb) then

            call nvtqvv_b1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,taut,sigma,tolnce,tstep,
     x        vircom,vircon)

          else

            call nvtqvv_b2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,taut,sigma,tolnce,tstep,
     x        vircom,vircon)
            
          endif

        elseif(keyens.eq.3) then
          
          if(.not.lcnb) then 

            call nvtqvv_h1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntshl,keyshl,chit,consv,conint,engke,
     x        engrot,taut,sigma,tolnce,tstep,vircom,vircon,chit_shl,
     x        sigma_shl)

          else

            call nvtqvv_h2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntshl,keyshl,chit,consv,conint,engke,
     x        engrot,taut,sigma,tolnce,tstep,vircom,vircon,chit_shl,
     x        sigma_shl)

          endif
          
        elseif(keyens.eq.4) then

          if(.not.lcnb) then

            call nptqvv_b1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,engke,engrot,press,taut,taup,sigma,
     x        tolnce,tstep,vircom,vircon,elrc,virlrc,virtot,volm)
            
          else

            call nptqvv_b2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,engke,engrot,press,taut,taup,sigma,
     x        tolnce,tstep,vircom,vircon,elrc,virlrc,virtot,volm)

          endif

        elseif(keyens.eq.5) then
          
          if(.not.lcnb) then 

            call nptqvv_h1
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,ntshl,keyshl,tstep,taut,taup,sigma,
     x        temp,chip,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircon,virtot,virlrc,vircom,volm,press,chit_shl,
     x        sigma_shl)

          else
      
            call nptqvv_h2
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,ntshl,keyshl,tstep,taut,taup,sigma,
     x        temp,chip,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircom,vircon,virtot,virlrc,volm,press,chit_shl,
     x        sigma_shl)
      
          endif
          
        elseif(keyens.eq.6) then

          if(.not.lcnb) then

            call nstqvv_b1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,engke,engrot,press,taut,taup,
     x        sigma,tolnce,tstep,vircom,vircon,elrc,virlrc,volm)

          else

            call nstqvv_b2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,engke,engrot,press,taut,taup,
     x        sigma,tolnce,tstep,vircom,vircon,elrc,virlrc,volm)

          endif

        elseif(keyens.eq.7) then

          if(.not.lcnb) then

            call nstqvv_h1
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,mode,ntshl,keyshl,tstep,taut,taup,
     x        sigma,temp,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircon,virlrc,vircom,volm,press,chit_shl,sigma_shl)

          else

            call nstqvv_h2
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,mode,ntshl,keyshl,tstep,taut,taup,
     x        sigma,temp,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircom,vircon,virlrc,volm,press,chit_shl,sigma_shl)
            
          endif

        else

c     invalid option

          call error(idnode,430)

        endif

      endif

c     check on convergence of pmf-shake

      if(ntpmf.gt.0) then

        if(mxnode.gt.1) call gstate(safep)
        if(.not.safep) call error(idnode,438)

      endif    

c     check on convergence of shake

      if(ntcons.gt.0) then

        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe) call error(idnode,105)

      endif    

c     eliminate "flying ice cube" in long simulations (Berendsen)
      
      if(lnfic.and.(keyens.eq.2.or.keyens.eq.4.or.keyens.eq.6))then
        
        if(mod(nstep,nofic).eq.0)then
          
          call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
          
        endif
        
      endif
      
      return
      end subroutine vv_integrate

      end module integrator_module
