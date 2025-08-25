let
    points = Observable(Point2f[])
    cond = Threads.Condition()
    global function selectbounds(fig::Figure, ax, x=[], y=[], z=[]; init=Point2f[], markersize=18, markercolor=:red)
        GLMakie.activate!()

        resize!(points[], length(init))
        points[] = init

        lines!(ax, points)
        p = scatter!(ax, points; markersize=markersize, color=markercolor)

        dragging = false
        draggingidx = 0

        isheatmap = !isempty(x) && !isempty(y) && !isempty(z)
        if isheatmap
            # linecutfig = Figure()
            linecutax = Axis(fig[2, :], xlabel="x", ylabel="y")
            linecutp = scatter!(linecutax, [0], [0])
            highlightp = vspan!(linecutax, -1, 1; color=(:gray, 0.4))
            pickingx = Ref(0.0)
            dy = mean([y[i+1] - y[i] for i in eachindex(y)[1:end-1]])
        end

        on(events(fig).mousebutton, priority=2) do event
            if isheatmap && event.button == Mouse.right && event.action == Mouse.press
                plt, i = pick(fig)
                if plt === p
                    empty!(linecutax)
                    idx = argmin(abs.(x .- points[][i][1]))
                    pickingx[] = x[idx]
                    linecutp = scatter!(linecutax, y, z[:, idx]; markersize=1.2markersize, label="x = $(pickingx[])")
                    highlightp = vspan!(linecutax, points[][i][2] - dy / 2, points[][i][2] + dy / 2; color=(:gray, 0.4))
                    axislegend(linecutax)
                end
            elseif event.button == Mouse.left && event.action == Mouse.press
                plt, i = pick(fig)
                if Keyboard.d in events(fig).keyboardstate
                    if plt == p
                        deleteat!(points[], i)
                        notify(points)
                        return Consume()
                    end
                elseif Keyboard.a in events(fig).keyboardstate
                    push!(points[], mouseposition(ax))
                    sort!(points[])
                    notify(points)
                    return Consume()
                elseif isheatmap && plt === linecutp
                    xidx = findfirst(p -> p[1] ≈ pickingx[], points[])
                    # @info xidx, pickingx[], points[]
                    if isnothing(xidx)
                        push!(points[], Point2f(pickingx[], y[i]))
                        sort!(points[])
                    else
                        points[][xidx] = Point2f(pickingx[], y[i])
                    end
                    delete!(linecutax, highlightp)
                    highlightp = vspan!(linecutax, y[i] - dy / 2, y[i] + dy / 2; color=(:gray, 0.4))
                    notify(points)
                    return Consume()
                else
                    dragging = plt == p
                    draggingidx = i
                    return Consume(dragging)
                end
            elseif event.action == Mouse.release
                dragging = false
                return Consume(false)
            end
            return Consume(false)
        end

        on(events(fig).mouseposition, priority=2) do mp
            if dragging
                points[][draggingidx] = mouseposition(ax)
                notify(points)
                return Consume(true)
            end
            return Consume(false)
        end

        display(fig)
        # on(events(fig).keyboardbutton) do event
        #     if event.action == Keyboard.press && event.key == Keyboard.enter
        #         lock(() -> Threads.notify(cond), cond)
        #     end
        # end
        # if isheatmap
        #     on(events(fig).mousebutton, priority=1) do event
        #         if event.button == Mouse.left && event.action == Mouse.press
        #             plt, i = pick(linecutfig)

        #         end
        #     end
        # end

        on(events(fig).window_open) do isopen
            isopen[] || lock(() -> Threads.notify(cond), cond)
        end

        lock(() -> wait(cond), cond)
        delete!(ax, p)

        # for p in points[]
        #     println((collect(p)...,))
        # end
        # @info "done"
        sleep(0.1)
        return copy(points[])
    end
end

# selectbounds_async(fig::Figure, ax; init=Point2f[]) = @async selectbounds(fig, ax; init=init)

# fig, ax, p = heatmap(rand(100, 100))
# @async selectbounds(fig, ax)