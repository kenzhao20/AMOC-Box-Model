#testing out a heaviside function

def Johnson2007(tau,A_GM,kappa,E,r,D0,St0,Sd0,Ss0,Sn0,Td0,a):
    
    #Function that:
        #Uses leapfrog method with an asselin time-filter to solve the equations (found in Helen Johnson et al.'s 2007 paper) governing the
        #circulation of salt and heat in a 4-box model of the AMOC. Solves for the salinity and temperature of the boxes and thermocline depth
    
        #The funciton then compresses the data-points and saves as a .npy file
    
    #Inputs:
        #Important parameters:
            #tau ~ O(0.1)      #Zonal wind stress [Nm-2]
            #A_GM ~ O(10^3)    #Eddy diffusion coefficient [m2s-1]
            #kappa ~ O(10^-5)  #Diapycnal mixing coefficient [m2s-1]
            #E ~ O(10^6)       #Freshwater evaporative flux [m3s-1]
            #r ~ O(10^7)       #Constant volume flux between thermohaline box and northern box from lateral mixing from ekman cells and wind-driven gyres [m3s-1]
        #Initial conditions:
            #Td0 ~ 5           #Initial deep ocean temperature [C]
            #Sd0 ~ 35          #Initial deep ocean salinity [psu]
            #St0 ~ 35          #Initial thermocline salinity [psu]
            #Ss0 ~ 35          #Initial southern salinity [psu]
            #Sn0 ~ 35          #Initial northern salinity [psu]
            #D0 ~ 5            #Initial thermocline depth [m]
        #Misc:
            #a ~ 0.1           #Asselin time filter parameter
    

    import os
    import numpy as np
    
    #Set time step

    dt=3e4 #~1/1000 years
    ti=0
    tf=2e11 #~6000 yrs

    t=np.arange(ti,tf,dt)

    sizet = len(t)
    
    #Create folder to hold data
    
    path = 'tau:'+str(tau)+' A_GM:'+str(A_GM)+' kappa:'+str(kappa)+' E:'+str(E)+' r:'+str(r)+' Td0:'+str(Td0)+' St0:'+str(St0)+' Sd0:'+str(Sd0)+' Ss0:'+str(Ss0)+' Sn0:'+str(Sn0)+' D0:'+str(D0)+' a:'+str(a)
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    #Set parameters not under investigation

    A = 2.6e14 # Horizontal area of upwelling (of thermohaline box) [m2]
    Lx = 3e7 # zonal extent of southern ocean (length?) [m]
    fS = -1e-4 # Coriolis parameter in southern ocean box [s-1]
    Ly = 1e6 # meridional extent of southern ocean frontal region (width?) [m]
    g = 9.8 #acceleration due to gravity [ms-2]
    fN = 1e-4 #Coriolis parameter in northern ocean box [s-1]

    Vn = 3e15 #Volume of northern ocean boxes [m3]
    Vs = 9e15 #Volume of southern ocean box [m3]
    Vtot = 1.2e18 #Volume of total ocean [m3]

    #Temperature of northern, southern, and thermocline boxes held fixed.
    Tn = 5 #Temperature of nothern ocean[C]
    Ts = 5 #Temperature of southern ocean[C]
    Tt = 25 #Temperature of thermohaline ocean[C]

    S0 = 35 #reference salinity [psu]
    T0 = 5 #Reference temperature [C]
    rho0 = 1027.5 # Reference density [kgm-3]

    alpha = 0.0002 #coefficient of fractional density change per change in temperature from reference T[C-1]
    beta = 0.0008 #coefficient of fractional density change per change in salinity from reference S[psu-1]
    
    
    #Functions to calculate said variables
    def rho(rho0,alpha,beta,T,T0,S,S0):  
        #Calculates density from temperature and salinity
        density = rho0*(1-alpha*(T-T0)+beta*(S-S0))
        return density

    def Vthermo(A,D):
        #Calculates volume of thermocline box
        Vt = A*D
        return Vt

    def Vdeep(Vtot,Vn,Vs,Vt):
        #Calculates volume of the deep ocean box
        Vd=Vtot-Vn-Vs-Vt
        return Vd

    def q_Ekman(tau,Lx,rho0,fS): 
        #Parameterises ekman flux from southern -> thermohaline
        qEk = tau*Lx/rho0/abs(fS)
        return qEk

    def q_Eddy(A_GM,D,Lx,Ly):
        #Parameterises eddy fluxes from thermohaline -> southern
        qEd = A_GM*D*Lx/Ly
        return qEd

    def q_South(qEk,qEd):
        #Net flux from southern -> thermohaline
        qS = qEk-qEd
        return qS

    def q_Up(kappa,A,D):
        #Parameterises upwelling flux from deep -> thermohaline
        qU = kappa*A/D
        return qU

    def q_North(g,rhod,rhot,rho0,rhon,D,fN):
        
        m = np.heaviside(rhon-rhod,1)
        
        #m=0.5*(1+np.tanh(k*(rhon-rhod)))
        
        qN = m*g*(rhot-rhod)*D**2/(2*fN*rho0)

        return qN

    #Create arrays of variables

    D = np.zeros(sizet)#Thermohaline depth [m]
    St = np.zeros(sizet) #Salinity of thermohaline box [psu]
    Sd = np.zeros(sizet) #Salinity in deep oceans [psu]
    Ss = np.zeros(sizet) #Salinity in southern oceans [psu]
    Sn = np.zeros(sizet) #Salinity in northern oceans [psu]
    Td = np.zeros(sizet) #Deep ocean temperature

    #Set initial conditions

    D[0] = D0
    
    St[0] = St0
    Sd[0] = Sd0
    Ss[0] = Ss0
    Sn[0] = Sn0

    Td[0] = Td0

    #Use Explicit Finite difference to generate the 1st and 2nd time-indexed slot:

    for i in range(2):
        #Find densities

        rhot = rho(rho0,alpha,beta,Tt,T0,St[i],S0)
        rhod = rho(rho0,alpha,beta,Td[i],T0,Sd[i],S0)
        rhos = rho(rho0,alpha,beta,Ts,T0,Ss[i],S0)
        rhon = rho(rho0,alpha,beta,Tn,T0,Sn[i],S0)

        #Find volume
        Vt = Vthermo(A,D[i])
        Vd = Vdeep(Vtot,Vn,Vs,Vt)

        #Find fluxes
        qEk = q_Ekman(tau,Lx,rho0,fS)
        qEd = q_Eddy(A_GM,D[i],Lx,Ly)
        qS = q_South(qEk,qEd)
        qU = q_Up(kappa,A,D[i])
        qN = q_North(g,rhod,rhot,rho0,rhon,D[i],fN)

        Vflux = qEk-qEd+qU-qN #Volume flux into the thermohaline box2

        #Use if, then for now. Will upgrade to a less computationally expensive method later

        #m = 0.5*(1+np.tanh(k*(qS)))
        
        m = np.heaviside(qS,1)

        D[i+1] = D[i] + dt * Vflux/A
        St[i+1] = St[i] + dt * (qU*Sd[i]+2*E*S0-qN*St[i]+(qS*(m*Ss[i]+(1-m)*St[i]))+r*(Sn[i]-St[i])-St[i]*Vflux)/Vt
        Sn[i+1] = Sn[i] + dt * (qN*(St[i]-Sn[i])-E*S0+r*(St[i]-Sn[i]))/Vn
        Sd[i+1] = Sd[i] + dt * (qN*Sn[i]-qU*Sd[i]-qS*(m*Sd[i]+(1-m)*Ss[i])+Sd[i]*Vflux)/Vd 
        Ss[i+1] = Ss[i] + dt * (qS*(m*(Sd[i]-Ss[i])+(1-m)*(Ss[i]-St[i]))-E*S0)/Vs
        Td[i+1] = Td[i] + dt * (qN*Tn-qU*Td[i]-(qS*(m*Td[i]+(1-m)*Ts))+Td[i]*Vflux)/Vd

    #Now use leapfrog or centred difference scheme and Asselin time filter to generate all other entries

    for i in range(3,sizet):
        #Find densities

        rhot = rho(rho0,alpha,beta,Tt,T0,St[i-1],S0)
        rhod = rho(rho0,alpha,beta,Td[i-1],T0,Sd[i-1],S0)
        rhos = rho(rho0,alpha,beta,Ts,T0,Ss[i-1],S0)
        rhon = rho(rho0,alpha,beta,Tn,T0,Sn[i-1],S0)

        #Find volume
        Vt = Vthermo(A,D[i-1])
        Vd = Vdeep(Vtot,Vn,Vs,Vt)

        #Find fluxes
        qEk = q_Ekman(tau,Lx,rho0,fS)
        qEd = q_Eddy(A_GM,D[i-1],Lx,Ly)
        qS = q_South(qEk,qEd)
        qU = q_Up(kappa,A,D[i-1])
        qN = q_North(g,rhod,rhot,rho0,rhon,D[i-1],fN)

        Vflux = qEk-qEd+qU-qN #Volume flux into the thermohaline box2

        #Use if, then for now. Will upgrade to a less computationally expensive method later
        
        m = 0.5*(1+np.tanh(k*(qS)))
        

        D[i] = D[i-2] + a*(D[i-3]-2*D[i-2]+D[i-1]) + 2 * dt/A * Vflux
        St[i] = St[i-2] + a*(St[i-3]-2*St[i-2]+St[i-1]) + 2 * dt/Vt * (qU*Sd[i-1]+2*E*S0-qN*St[i-1]+(qS*(m*Ss[i-1]+(1-m)*St[i-1]))+r*(Sn[i-1]-St[i-1])-St[i-1]*Vflux)
        Sn[i] = Sn[i-2] + a*(Sn[i-3]-2*Sn[i-2]+Sn[i-1]) + 2 * dt/Vn * (qN*(St[i-1]-Sn[i-1])-E*S0+r*(St[i-1]-Sn[i-1]))
        Sd[i] = Sd[i-2] + a*(Sd[i-3]-2*Sd[i-2]+Sd[i-1]) + 2 * dt/Vd * (qN*Sn[i-1]-qU*Sd[i-1]-qS*(m*Sd[i-1]+(1-m)*Ss[i-1])+Sd[i-1]*Vflux)
        Ss[i] = Ss[i-2] + a*(Ss[i-3]-2*Ss[i-2]+Ss[i-1]) + 2 * dt/Vs * (qS*(m*(Sd[i-1]-Ss[i-1])+(1-m)*(Ss[i-1]-St[i-1]))-E*S0)
        Td[i] = Td[i-2] + a*(Td[i-3]-2*Td[i-2]+Td[i-1]) + 2 * dt/Vd * (qN*Tn-qU*Td[i-1]-(qS*(m*Td[i-1]+(1-m)*Ts))+Td[i-1]*Vflux)


    #Calculate the densities at each time from the arrays:

    rhot = rho(rho0,alpha,beta,Tt,T0,St,S0)
    rhod = rho(rho0,alpha,beta,Td,T0,Sd,S0)
    rhos = rho(rho0,alpha,beta,Ts,T0,Ss,S0)
    rhon = rho(rho0,alpha,beta,Tn,T0,Sn,S0)

    #Calculate fluxes at each time from the arrays:
    qEk = q_Ekman(tau,Lx,rho0,fS)
    qEd = q_Eddy(A_GM,D,Lx,Ly)
    qS = q_South(qEk,qEd)
    qU = q_Up(kappa,A,D)
    #qN = q_North(g,rhod,rhot,rho0,rhon,D,fN)


    #Save parameters used
    parameters = np.array([('tau',tau),('A_GM',A_GM),('kappa',kappa),('E',E),('r',r),('asselin time filter',a)])

    np.save(path+'/parameters',parameters)

    #Compress data files and save
    ndp = 400 #Number of data points
    entries = round(sizet/ndp)

    np.save(path+'/t',t[::entries])

    np.save(path+'/D',D[::entries])

    np.save(path+'/St',St[::entries])
    np.save(path+'/Sn',Sn[::entries])
    np.save(path+'/Sd',Sd[::entries])
    np.save(path+'/Ss',Ss[::entries])

    np.save(path+'/Td',Td[::entries])

    np.save(path+'/rhot',rhot[::entries])
    np.save(path+'/rhod',rhod[::entries])
    np.save(path+'/rhos',rhos[::entries])
    np.save(path+'/rhon',rhon[::entries])

    np.save(path+'/qEk',qEk)
    np.save(path+'/qEd',qEd[::entries])
    np.save(path+'/qS',qS[::entries])
    np.save(path+'/qU',qU[::entries])