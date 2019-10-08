using PGFPlots
using Reel
using Colors

# Setup plotting for animinations
# Set default axis and label styles
Reel.extension(m::MIME"image/svg+xml") = "svg"
Reel.set_output_type("gif") # may be necessary for use in IJulia


function plotReachable(sets,extent,time;zmax=360,timeOffset=0,cmCustom=true,desiredKey=nothing,psi=nothing,nbin=150,width="18cm",height="18cm",title=-1,xlabel="Downrange (kft)",ylabel="Crossrange (kft)",filename="Pics/HCAS_ReachablePlot.png",savePlot=false,saveName=saveName="Pics/HCAS_ReachablePlot.tex",isRow=false,isLast=false,keepReg=nothing)
    #= Plot reachable set
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        time (int): Time (tau) value to plot (Should be a key in sets)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom Jet colormap
        desiredKey (string): Previous advisory sequence to plot (should be a key of sets[time]).
            If nothing, the reachable sets of all previous advisory sequences are combined
        psi (float): Show reachable cells that contain a given psi value.
            If nothing, all reachable cells are used
        nbin (int): Number of bins for plotting resolution
        width (string): Width of plot
        height (string): height of plot
        title (string): Title of plot. If title==-1, use a default title
        xlabel (string): X-axis label
        ylabel (string): Y-axis label
        filename (string): Filename for saving png file (if saving)
        savePlot (bool): True if saving the plot as .tex and .png files
        saveName (string): Name for the .tex file
        isRow (bool): True if this plot is a part of a larger row of plots
        isLast (bool): True if this is the last plot in an animation
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        g (PGFPlots.Axis): Plot of reachable regions
    =#
    set = BitArray(undef,NUMREGIONS)
    if desiredKey != nothing
        set .= sets[time-timeOffset][desiredKey]
    else
        set .= false
        for (k,v) = sets[time-timeOffset]
            set .|= v
        end
    end
    
    if !cmCustom
        zmax=1
    end
    
    function getXY(x,y)
        x*=1000
        y*=1000
        
        
        if y>YS[end-1]
            y=YS[end-1]
        elseif y<YS[1]
            y=YS[1]
        end
        if x > XS[end-1]
            x=XS[end-1]
        elseif x<XS[1]
            x=XS[1]
        end
        
        xInd = findall(XS.>x)[1]-1
        yInd = findall(YS.>y)[1]-1
        sumFound=0
        if psi !=nothing
            if psi>PSIS[end-1]
                psi=PSIS[end-1]
            elseif psi<PSIS[1]
                psi = PSIS[1]
            end
            pInd = findall(PSIS.>psi)[1]-1
            ind = indicesToIndex(pInd,yInd,xInd)
            if (keepReg != nothing) && (!keepReg[ind])
                return zmax
            end
            sumFound = set[ind]
        else
            inds = indicesToIndex_Vector(1:NUMP,yInd*ones(Int32,NUMP),xInd*ones(Int32,NUMP))
            if (keepReg != nothing) && (!keepReg[inds[1]])
                return zmax
            end
            sumFound = sum([set[i] for i in inds])
        end
        return sumFound
    end
    
    if !savePlot && !isRow
        filename=nothing
    end
    
    if !isRow
        PGFPlots.resetPGFPlotsPreamble()
        PGFPlots.pushPGFPlotsPreamble("\\pgfplotsset{every axis/.append style={title style={font=\\huge},label style={font=\\huge},tick label style={font=\\Large}}}")
    end
        
    xmin, xmax, ymin, ymax = extent
    if title==-1
        title=@sprintf("Time to Coaltitude: %d s",max(0,time))
        if time<=0
            title=@sprintf("Time at Coaltitude: %d s",-time)
        end
        if isLast
            title="Converged to Steady State"
        end
    end
    style=""
    if ylabel==""
        style = "yticklabels={,,}, scaled y ticks=false"
    end
    cm = ColorMaps.Named("Blues")
    if cmCustom
        cmColors = copy(ColorMaps.Named("Jet").colors)
        cmColors[1] = RGB(1.0,1.0,1.0)
        cm=ColorMaps.RGBArrayMap(cmColors)
    end
    g = Axis([
            Plots.Image(getXY,(xmin,xmax),(ymin,ymax),zmin=0,zmax=zmax,
                xbins=nbin,ybins=nbin,colormap=cm,colorbar=false,filename=filename),
            Plots.Circle(0,0,0.5,style="black,fill=red"),
        ],
        width=width,height=height,xlabel=xlabel,ylabel=ylabel,title=title,style=style)
    if savePlot
        save(saveName,g,include_preamble=false)
    else
        return g
    end
end


function plotReachableRow(sets,extent,times;nbin=200,zmax=360,timeOffset=0,cmCustom=true,savePlot=false,saveName="Pics/HCAS_ReachableRow.tex",keepReg=nothing)
    #= Plot a row of reachable plots, which are the reachable regions at different times
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        times (array of int): Times (taus) at which we want to plot (Should be keys in sets)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom Jet colormap
        nbin (int): Number of bins for plotting resolutiones
        saveName (string): Name for the .tex file
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        g (PGFPlots.GroupPlot): Row of reachable region plots
    =#
    PGFPlots.resetPGFPlotsPreamble()
    g = GroupPlot(length(times),1, groupStyle = "horizontal sep=0.35cm", style="height=4.4cm, width=4.4cm")
    width=nothing; height=nothing;ylabel="Crossrange (kft)";xlabel="Downrange (kft)"
    for time in times
        fn = nothing
        if savePlot
            fn = @sprintf("Pics/HCAS_ReachRow_%d.png",time)
        end
        push!(g,plotReachable(sets,extent,time,nbin=nbin,width=width,height=height,title=-1,xlabel=xlabel,ylabel=ylabel,filename=fn,zmax=zmax,timeOffset=timeOffset,cmCustom=cmCustom,isRow=true,keepReg=keepReg))
        ylabel=""
    end
    if savePlot
        save(saveName,g,include_preamble=false)
        return
    else
        return g
    end
end

function animateReach(sets,extent,timeMax,timeMin,save_name;timeOffset=0,cmCustom=true,zmax=360,nbin=150,keepReg=nothing)
    #= Create an animation of reachable set plots. Begins with timeMax and counts down to timeMin
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        timeMin (int): Minimum time (tau) value to plot (Should be a key in  sets)
        timeMax (int): Maximum time (tau) value to plot (Should be a key in  sets)
        save_name (string): Name of gif file (without the .gif)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom Jet colormap
        nbin (int): Number of bins for plotting resolution
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        frames (Reel.Frames): Animation of reachable regions changing over time
    =#
    
    currentFrac = 0
    deltaFrac = 0.1
    frames = Frames(MIME("image/svg+xml"), fps=10)
    for frame in timeMax:-1:timeMin
        if (timeMax-frame)/(timeMax-timeMin)>=currentFrac
            println("Completed "*string(round(currentFrac*100))*"%")
            currentFrac+=deltaFrac
        end
        push!(frames,plotReachable(sets,extent,frame,timeOffset=timeOffset,cmCustom=cmCustom,zmax=zmax,nbin=nbin,isLast=frame==timeMin-1,keepReg=keepReg))
    end
    write(save_name*".gif",frames)
    frames
end