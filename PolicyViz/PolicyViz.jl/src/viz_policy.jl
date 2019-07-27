export viz_policy, viz_policy_full, viz_values, percCorrect

"""
Helper function to draw aircraft
"""
function getACString(theta,x,y,fill,draw,width=2.0)
    return "\\node[aircraft top,fill="*fill*",draw="*draw*", minimum width="*string(width)*"cm,rotate="*string(theta)*",scale = 0.35] at (axis cs:"*string(x)*", "*string(y)*") {};"
end

function percCorrect(nnetPath,tablePath)
    

    nnet = read_network(nnetPath);
    Q = h5open(tablePath, "r") do file
        read(file, "y")
    end
    X = h5open(tablePath, "r") do file
        read(file, "X")
    end

    rnges = h5open(tablePath, "r") do file
        read(file, "ranges")
    end
    means = h5open(tablePath, "r") do file
        read(file, "means")
    end

    Q = Q.*rnges[end].+means[end]
    X = X.*nnet.ranges[1:end-1].+nnet.means[1:end-1]
    nnet_policy = argmax(evaluate_network_multiple(nnet,X),dims=[1])
    table_policy = argmax(Q,dims=[1])
    return sum(nnet_policy.==table_policy)/length(nnet_policy)
end

function viz_values(;nnetPath::AbstractString="",tablePath::AbstractString="",batch_size::Int=-1)
    if nnetPath=="" && tablePath==""
        println("Please provide at least either a path to the neural network file or training data file")
        return
    end
    
    # Use aircraft shapes tikz style package
    pushPGFPlotsPreamble("\\usepackage{aircraftshapes}")
    
    # Load network if given
    nnet = nothing
    if nnetPath!=""
        nnet = read_network(nnetPath);
    end
    
    # Load table if given
    Q = nothing
    if tablePath!=""
        Q = h5open(tablePath, "r") do file
            read(file, "y")
        end

        rnges = h5open(tablePath, "r") do file
            read(file, "ranges")
        end
        means = h5open(tablePath, "r") do file
            read(file, "means")
        end

        Q = Q.*rnges[end].+means[end]
    end

    grid  = RectangleGrid(RANGES,THETAS,PSIS,OWNSPEEDS,INTRSPEEDS)

    # Create manipulatable plotting tool
    @manipulate for psi  = convert(Array{Int32,1},round.(rad2deg.(PSIS))),
        ownspeed = OWNSPEEDS,
        intrspeed = INTRSPEEDS,
        zoom = 1.0,
        nbin = [100,150,200,250,1000],
        add_cocPenalty = [false,true],
        add_online = [false,true],
        savePlot = [false,true],
        xshift = 0.0,
        yshift = 0.0,
        xscale = 1.0,
        yscale = 1.0,
        ra     = [0,1,2,3,4],
        zmin   = -10.0,
        zmax   = 0.0
        
        # Ensure that zoom and scale factors don't result in dividing by zero or negative numbers
        zoom = zoom <0.1 ? 0.1 : zoom
        xscale = xscale <0.1 ? 0.1 : xscale
        yscale = yscale <0.1 ? 0.1 : yscale
        
        # Compute neural network Q-values if a neural network is given
        if nnet != nothing
            
            #Load table with the inputs needed to plot the heat map
            inputsNet= zeros(5,nbin*nbin)    
            ind = 1
            for i=range((-1*RANGEMAX/zoom/xscale- xshift)*1000,stop=(RANGEMAX/zoom/xscale- xshift)*1000,length=nbin)
                for j=range((-1*RANGEMAX/zoom/yscale - yshift)*1000,stop=(RANGEMAX/zoom/yscale- yshift)*1000,length=nbin)
                    r = sqrt(i^2+j^2)
                    th = atan(j,i)
                    inputsNet[:,ind] = [r,th,deg2rad.(psi),ownspeed,intrspeed]
                    ind = ind+1
                end
            end

            #Calculate all of the Q values from the input array
            q_nnet = zeros(5,nbin*nbin);
            ind = 1
            batch_size = batch_size<0 ? nbin*nbin : batch_size

            # Compute network outputs in batches
            while ind+batch_size<nbin*nbin
                input = inputsNet[:,ind:(ind+batch_size-1)]
                q_nnet[:,ind:(ind+batch_size-1)] = evaluate_network_multiple(nnet,input)
                ind=ind+batch_size
            end
            input = inputsNet[:,ind:end]
            q_nnet[:,ind:end] = evaluate_network_multiple(nnet,input)
        end

        
        # Q Table Heat Map
        function get_heat1(x::Float64, y::Float64)
            r = sqrt(x^2+y^2)*1000
            th = atan(y,x)
            bel = get_belief([r,th,deg2rad.(psi),ownspeed,intrspeed],grid,false)
            qvals = Q[:,bel.rowval[1]]
            
            return qvals[ra+1]
        end # function get_heat
        
        
       ind = 1       
       #Neural Net Heat Map
       function get_heat2(x::Float64, y::Float64)  
            qvals = q_nnet[:,ind]
            ind += 1
            r = sqrt(x^2+y^2)*1000
            th = atan(y,x)
            
            return qvals[ra+1]
       end # function get_heat2

        numPlots = 2
        numPlots = nnet==nothing ? numPlots-1 : numPlots
        numPlots = Q==nothing ? numPlots-1 : numPlots
        
        g = GroupPlot(numPlots, 1, groupStyle = "horizontal sep=3.5cm")
        
        # Plot Q-table, if possible
        if Q != nothing
            push!(g, Axis([
                        Plots.Image(get_heat1, (-1*RANGEMAX/zoom/xscale - xshift, RANGEMAX/zoom/xscale- xshift), 
                            (-1*RANGEMAX/zoom/yscale- yshift, RANGEMAX/zoom/yscale- yshift), 
                            zmin = zmin, zmax = zmax,
                            xbins = nbin, ybins = nbin,
                            colormap = ColorMaps.Named("Jet"), colorbar=true),
                        Plots.Command(getACString(0.0,0.0,0.0,"black","white")),
                        Plots.Command(getACString(psi,RANGEMAX*0.85/zoom/xscale - xshift,RANGEMAX*0.85/zoom/yscale- yshift,"red","white"))
                        ], width="8cm", height="8cm", xlabel="Downrange (kft)", ylabel="Crossrange (kft)", title="Training Data"))
        end
        
        # Plot neural network, if possible
        if nnet != nothing
            push!(g, Axis([
                        Plots.Image(get_heat2, (-1*RANGEMAX/zoom/xscale - xshift, RANGEMAX/zoom/xscale- xshift), 
                            (-1*RANGEMAX/zoom/yscale- yshift, RANGEMAX/zoom/yscale- yshift), 
                            zmin = zmin, zmax = zmax,
                            xbins = nbin, ybins = nbin,
                            colormap = ColorMaps.Named("Jet"), colorbar=true),
                        Plots.Command(getACString(0.0,0.0,0.0,"black","white")),
                        Plots.Command(getACString(psi,RANGEMAX*0.85/zoom/xscale - xshift,RANGEMAX*0.85/zoom/yscale - yshift,"red","white"))
                        ], width="8cm", height="8cm", xlabel="Downrange (kft)", ylabel="Crossrange (kft)", title="Neural Network"))
        end
        if savePlot
            PGFPlots.save("ValuePlot.tex", g, include_preamble=true)
        else           
            return g
        end
        return nothing
    end # for p_int, v0, v1, pa, ta, etc
end # function viz_policy
    
"""
Visualize the neural network and/or table policies
"""
function viz_policy(;nnetPath::AbstractString="",tablePath::AbstractString="",batch_size::Int=-1,nnetRect=false)
    if nnetPath=="" && tablePath==""
        println("Please provide at least either a path to the neural network file or training data file")
        return
    end
    
    # Use aircraft shapes tikz style package
    resetPGFPlotsPreamble()
    pushPGFPlotsPreamble("\\usepackage{aircraftshapes}")
    
    # Load network if given
    nnet = nothing
    if nnetPath!=""
        nnet = read_network(nnetPath);
    end
    
    # Load table if given
    Q = nothing
    if tablePath!=""
        Q = h5open(tablePath, "r") do file
            read(file, "y")
        end

        rnges = h5open(tablePath, "r") do file
            read(file, "ranges")
        end
        means = h5open(tablePath, "r") do file
            read(file, "means")
        end

        Q = Q.*rnges[end].+means[end]
    end

    grid  = RectangleGrid(RANGES,THETAS,PSIS,OWNSPEEDS,INTRSPEEDS)
    
    COC = RGB(1.,1.,1.) # white
    SR = RGB(.0,.0,.5) # navy
    SL = RGB(.0,.600,.0) # green
    WR = RGB(.5,.5,.5) # grey
    WL = RGB(.7,.9,.0) # neon green
    ra_colors = [SL,WL,COC,WR,SR]
    bg_colors = [COC]
    
    # Create scatter plot classes for color key
    sc_string = "{"
    for i=0:4
        define_color("ra_$i",  ra_colors[i+1])
        if i==2
            sc_string *= "ra_$i={mark=square, style={black, mark options={fill=ra_$i}, mark size=6}},"
        else
            sc_string *= "ra_$i={style={ra_$i, mark size=6}},"
        end
    end
    
    # Color key as a scatter plot
    sc_string=sc_string[1:end-1]*"}"
    xx = [-1.5,-1.5,-1.5, -1.5, -1.5]
    yy = [1.65,1.15,0.65, 0.15, -0.35]
    zz = ["ra_0","ra_1","ra_2","ra_3","ra_4"]
    sc = string(sc_string)
    
    # Create manipulatable plotting tool
    @manipulate for psi  = convert(Array{Int32,1},round.(rad2deg.(PSIS))),
        ownspeed = OWNSPEEDS,
        intrspeed = INTRSPEEDS,
        zoom = 1.0,
        nbin = [100,150,200,250,1000],
        savePlot = [false,true],
        xshift = 0.0,
        yshift = 0.0,
        xscale = 1.0,
        yscale = 1.0
        
        # Ensure that zoom and scale factors don't result in dividing by zero or negative numbers
        zoom = zoom <0.1 ? 0.1 : zoom
        xscale = xscale <0.1 ? 0.1 : xscale
        yscale = yscale <0.1 ? 0.1 : yscale
        
        # Compute neural network Q-values if a neural network is given
        if nnet != nothing
            
            #Load table with the inputs needed to plot the heat map
            inputsNet= zeros(5,nbin*nbin)    
            ind = 1
            for i=range((-1*RANGEMAX/zoom/xscale- xshift)*1000,stop=(RANGEMAX/zoom/xscale- xshift)*1000,length=nbin)
                for j=range((-1*RANGEMAX/zoom/yscale - yshift)*1000,stop=(RANGEMAX/zoom/yscale- yshift)*1000,length=nbin)
                    r = sqrt(i^2+j^2)
                    th = atan(j,i)
                    if nnetRect
                        inputsNet[:,ind] = [i,j,deg2rad.(psi),ownspeed,intrspeed]
                    else
                        inputsNet[:,ind] = [r,th,deg2rad.(psi),ownspeed,intrspeed]
                    end
                    ind = ind+1
                end
            end

            #Calculate all of the Q values from the input array
            q_nnet = zeros(5,nbin*nbin);
            ind = 1
            batch_size = batch_size<0 ? nbin*nbin : batch_size

            # Compute network outputs in batches
            while ind+batch_size<nbin*nbin
                input = inputsNet[:,ind:(ind+batch_size-1)]
                q_nnet[:,ind:(ind+batch_size-1)] = evaluate_network_multiple(nnet,input)
                ind=ind+batch_size
            end
            input = inputsNet[:,ind:end]
            q_nnet[:,ind:end] = evaluate_network_multiple(nnet,input)
        end

        
        # Q Table Heat Map
        function get_heat1(x::Float64, y::Float64)
            r = sqrt(x^2+y^2)*1000
            if r>RANGEMAX*1000
                return 0.0
            end
            th = atan(y,x)
            bel = get_belief([r,th,deg2rad.(psi),ownspeed,intrspeed],grid,false)
            qvals = Q[:,bel.rowval[1]]
            return rad2deg.(ACTIONS[findmax(qvals)[2]])
        end # function get_heat
        
        
       ind = 1       
       #Neural Net Heat Map
       function get_heat2(x::Float64, y::Float64)  
            qvals2 = q_nnet[:,ind]
            ind += 1
            r = sqrt(x^2+y^2)*1000
            th = atan(y,x)
            if r>RANGEMAX*1000
                return 0.0
            end
            return rad2deg.(ACTIONS[findmax(qvals2)[2]])
       end # function get_heat2

        numPlots = 3
        numPlots = nnet==nothing ? numPlots-1 : numPlots
        numPlots = Q==nothing ? numPlots-1 : numPlots
        numPlots = savePlot ? numPlots-1 : numPlots
        
        g = GroupPlot(numPlots, 1, style="height={8cm}, width={8cm}",groupStyle = "horizontal sep=2cm")
        
        # Plot Q-table, if possible
        if Q != nothing
            push!(g, Axis([
                        Plots.Image(get_heat1, (-1*RANGEMAX/zoom/xscale - xshift, RANGEMAX/zoom/xscale- xshift), 
                            (-1*RANGEMAX/zoom/yscale- yshift, RANGEMAX/zoom/yscale- yshift), 
                            zmin = -3, zmax = 3,
                            xbins = nbin, ybins = nbin,
                            colormap = ColorMaps.RGBArrayMap(ra_colors), colorbar=false),
                        Plots.Command(getACString(0.0,0.0,0.0,"black","white")),
                        Plots.Command(getACString(psi,RANGEMAX*0.82/zoom/xscale - xshift,RANGEMAX*0.82/zoom/yscale- yshift,"red","black"))
                        ], xlabel="Downrange (kft)", ylabel="Crossrange (kft)", title="Training Data"))
        end
        
        # Plot neural network, if possible
        if nnet != nothing
            push!(g, Axis([
                        Plots.Image(get_heat2, (-1*RANGEMAX/zoom/xscale - xshift, RANGEMAX/zoom/xscale- xshift), 
                            (-1*RANGEMAX/zoom/yscale- yshift, RANGEMAX/zoom/yscale- yshift), 
                            zmin = -3, zmax = 3,
                            xbins = nbin, ybins = nbin,
                            colormap = ColorMaps.RGBArrayMap(ra_colors), colorbar=false),
                        Plots.Command(getACString(0.0,0.0,0.0,"black","white")),
                        Plots.Command(getACString(psi,RANGEMAX*0.82/zoom/xscale - xshift,RANGEMAX*0.82/zoom/yscale - yshift,"red","black"))
                        ], xlabel="Downrange (kft)", ylabel="Crossrange (kft)", title="Neural Network"))
        end
        
        # Create Color Key
        f = (x,y)->x # Dummy function for background white image
        push!(g, Axis([
            Plots.Image(f, (-2,2), (-2,2),colormap = ColorMaps.RGBArrayMap(bg_colors),colorbar=false),
            Plots.Scatter(xx, yy, zz, scatterClasses=sc),
            Plots.Node("SL ",0.25,0.915,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("WL ",0.25,0.790,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("COC",0.25,0.665,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("WR ",0.25,0.540,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("SR ",0.25,0.415,style="black,anchor=west", axis="axis description cs"),
            ],style="xshift=-1.4cm",hideAxis =true))
        
        # Save plot if desired. Don't return g if want to save, since returning g will cleanup and delete the png images.
        if savePlot
            PGFPlots.save("PolicyPlot.tex", g, include_preamble=true)
        else           
            return g
        end
        return nothing
    end # for p_int, v0, v1, pa, ta, etc
end # function viz_policy

"""
Visualize the neural network and/or table policies
"""
function viz_policy_full(;nnetPath::AbstractString="",useTable=false,batch_size::Int=-1)
    if nnetPath=="" && !useTable
        println("Please provide at least either a path to the neural network file or useTable")
        return
    end
    
    # Use aircraft shapes tikz style package
    pushPGFPlotsPreamble("\\usepackage{aircraftshapes}")
    
    # Load network if given
    nnet = nothing
    if nnetPath!=""
        nnet = read_network(nnetPath);
    end
    
    # Load table if given
    Q = nothing
    pra_current = -1
    tau_current = -1

    grid  = RectangleGrid(RANGES,THETAS,PSIS,OWNSPEEDS,INTRSPEEDS)
    
    COC = RGB(1.,1.,1.) # white
    SR = RGB(.0,.0,.5) # navy
    SL = RGB(.0,.600,.0) # green
    WR = RGB(.5,.5,.5) # grey
    WL = RGB(.7,.9,.0) # neon green
    ra_colors = [SL,WL,COC,WR,SR]
    bg_colors = [COC]
    
    # Create scatter plot classes for color key
    sc_string = "{"
    for i=0:4
        define_color("ra_$i",  ra_colors[i+1])
        if i==2
            sc_string *= "ra_$i={mark=square, style={black, mark options={fill=ra_$i}, mark size=6}},"
        else
            sc_string *= "ra_$i={style={ra_$i, mark size=6}},"
        end
    end
    
    # Color key as a scatter plot
    sc_string=sc_string[1:end-1]*"}"
    xx = [-1.5,-1.5,-1.5, -1.5, -1.5]
    yy = [1.65,1.15,0.65, 0.15, -0.35]
    zz = ["ra_0","ra_1","ra_2","ra_3","ra_4"]
    sc = string(sc_string)
    
    # Create manipulatable plotting tool
    @manipulate for psi  = convert(Array{Int32,1},round.(rad2deg.(PSIS))),
        ownspeed = OWNSPEEDS,
        intrspeed = INTRSPEEDS,
        tau = [0,5,10,15,20,30,40,60],
        pra = [0,1,2,3,4],
        zoom = 1.0,
        nbin = [100,150,200,250,1000],
        savePlot = [false,true],
        xshift = 0.0,
        yshift = 0.0,
        xscale = 1.0,
        yscale = 1.0
        
        # Ensure that zoom and scale factors don't result in dividing by zero or negative numbers
        zoom = zoom <0.1 ? 0.1 : zoom
        xscale = xscale <0.1 ? 0.1 : xscale
        yscale = yscale <0.1 ? 0.1 : yscale
        
        
        if (useTable) && ((pra_current!=pra) || (tau_current!=tau))
            tablePath = @sprintf("/raid/kjulian3/HCAS/TrainingData/HCAS_oneSpeed_TrainingData_v1_pra%d_tau%02d.h5",pra,tau)
            Q = h5open(tablePath, "r") do file
                read(file, "y")
            end

            rnges = h5open(tablePath, "r") do file
                read(file, "ranges")
            end
            means = h5open(tablePath, "r") do file
                read(file, "means")
            end

            Q = Q.*rnges[end].+means[end]
            pra_current = pra
            tau_current = tau
        end
        
        # Compute neural network Q-values if a neural network is given
        if nnet != nothing
            
            #Load table with the inputs needed to plot the heat map
            inputsNet= zeros(7,nbin*nbin)    
            ind = 1
            for i=range((-1*RANGEMAX/zoom/xscale- xshift)*1000,stop=(RANGEMAX/zoom/xscale- xshift)*1000,length=nbin)
                for j=range((-1*RANGEMAX/zoom/yscale - yshift)*1000,stop=(RANGEMAX/zoom/yscale- yshift)*1000,length=nbin)
                    r = sqrt(i^2+j^2)
                    th = atan(j,i)
                    inputsNet[:,ind] = [r,th,deg2rad.(psi),ownspeed,intrspeed,tau,pra]
                    ind = ind+1
                end
            end

            #Calculate all of the Q values from the input array
            q_nnet = zeros(5,nbin*nbin);
            ind = 1
            batch_size = batch_size<0 ? nbin*nbin : batch_size

            # Compute network outputs in batches
            while ind+batch_size<nbin*nbin
                input = inputsNet[:,ind:(ind+batch_size-1)]
                q_nnet[:,ind:(ind+batch_size-1)] = evaluate_network_multiple(nnet,input)
                ind=ind+batch_size
            end
            input = inputsNet[:,ind:end]
            q_nnet[:,ind:end] = evaluate_network_multiple(nnet,input)
        end

        
        # Q Table Heat Map
        function get_heat1(x::Float64, y::Float64)
            r = sqrt(x^2+y^2)*1000
            if r>RANGEMAX*1000
                return 0.0
            end
            th = atan(y,x)
            bel = get_belief([r,th,deg2rad.(psi),ownspeed,intrspeed],grid,false)
            qvals = Q[:,bel.rowval[1]]
            return rad2deg.(ACTIONS[findmax(qvals)[2]])
        end # function get_heat
        
        
       ind = 1       
       #Neural Net Heat Map
       function get_heat2(x::Float64, y::Float64)  
            qvals2 = q_nnet[:,ind]
            ind += 1
            r = sqrt(x^2+y^2)*1000
            th = atan(y,x)
            if r>RANGEMAX*1000
                return 0.0
            end
            return rad2deg.(ACTIONS[findmax(qvals2)[2]])
       end # function get_heat2

        numPlots = 3
        numPlots = nnet==nothing ? numPlots-1 : numPlots
        numPlots = Q==nothing ? numPlots-1 : numPlots
        numPlots = savePlot ? numPlots-1 : numPlots
        
        g = GroupPlot(numPlots, 1, groupStyle = "horizontal sep=2cm")
        
        # Plot Q-table, if possible
        if Q != nothing
            push!(g, Axis([
                        Plots.Image(get_heat1, (-1*RANGEMAX/zoom/xscale - xshift, RANGEMAX/zoom/xscale- xshift), 
                            (-1*RANGEMAX/zoom/yscale- yshift, RANGEMAX/zoom/yscale- yshift), 
                            zmin = -3, zmax = 3,
                            xbins = nbin, ybins = nbin,
                            colormap = ColorMaps.RGBArrayMap(ra_colors), colorbar=false),
                        Plots.Command(getACString(0.0,0.0,0.0,"black","white")),
                        Plots.Command(getACString(psi,RANGEMAX*0.85/zoom/xscale - xshift,RANGEMAX*0.85/zoom/yscale- yshift,"red","white"))
                        ], width="8cm", height="8cm", xlabel="Downrange (kft)", ylabel="Crossrange (kft)", title="Training Data"))
        end
        
        # Plot neural network, if possible
        if nnet != nothing
            push!(g, Axis([
                        Plots.Image(get_heat2, (-1*RANGEMAX/zoom/xscale - xshift, RANGEMAX/zoom/xscale- xshift), 
                            (-1*RANGEMAX/zoom/yscale- yshift, RANGEMAX/zoom/yscale- yshift), 
                            zmin = -3, zmax = 3,
                            xbins = nbin, ybins = nbin,
                            colormap = ColorMaps.RGBArrayMap(ra_colors), colorbar=false),
                        Plots.Command(getACString(0.0,0.0,0.0,"black","white")),
                        Plots.Command(getACString(psi,RANGEMAX*0.85/zoom/xscale - xshift,RANGEMAX*0.85/zoom/yscale - yshift,"red","white"))
                        ], width="8cm", height="8cm", xlabel="Downrange (kft)", ylabel="Crossrange (kft)", title="Neural Network"))
        end
        
        # Create Color Key
        f = (x,y)->x # Dummy function for background white image
        push!(g, Axis([
            Plots.Image(f, (-2,2), (-2,2),colormap = ColorMaps.RGBArrayMap(bg_colors),colorbar=false),
            Plots.Scatter(xx, yy, zz, scatterClasses=sc),
            Plots.Node("SL ",0.25,0.915,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("WL ",0.25,0.790,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("COC",0.25,0.665,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("WR ",0.25,0.540,style="black,anchor=west", axis="axis description cs"),
            Plots.Node("SR ",0.25,0.415,style="black,anchor=west", axis="axis description cs"),
            ],style="xshift=-1.4cm",width="5cm",height="8cm", hideAxis =true))
        
        # Save plot if desired. Don't return g if want to save, since returning g will cleanup and delete the png images.
        if savePlot
            PGFPlots.save("PolicyPlot.tex", g, include_preamble=true)
        else           
            return g
        end
        return nothing
    end # for p_int, v0, v1, pa, ta, etc
end # function viz_policy