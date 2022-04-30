# Elementy ze standardowej biblioteki
using Printf
using Plots

pyplot()
#plotlyjs()
#	====================================================================
#							Obliczenia
#   ====================================================================

#__precision = BigFloat
#setprecision(128)

__precision = Float64

eulerConstant = 0.57721566490153286060651209008240243104215933593992


# Konstruktor macierzy X*Y
function mat(dimX, dimY, initialV)
    mat1 = Array{__precision}(undef, dimX, dimY)
    
    for i in 1:dimX
        for j in 1:dimY
            mat1[i, j] = initialV
        end
    end
    
    return mat1
end


# n-ty wyraz ciągu gamma, sumowanie odbywa się od elementu najmniejszego tj. 1/n
function gamma(n)

    sum = zero(__precision)

    for i in 1:n
        sum = sum + 1/convert(__precision, (n - i + 1))
    end

    return  sum - log(convert(__precision, n))
end

#	====================================================================
#							Wyznaczanie wartości
#   ====================================================================

# algorytm regresji liniowej
function linear_regression(argumentsX, valuesY, verbose = false) # -> ax + b
    len = length(argumentsX)
	
	matX = mat(2, len, 1)
	matY = mat(1, len, 1)
	
	for i in 1:len
		matX[2, i] = argumentsX[i]
		matY[1, i] = valuesY[i]
	end

	a = (matX*transpose(matX))^-1* (matX) *transpose( matY)
	
    if verbose
        @printf("Przybliżenie a: %f\nb: %f\n", a[2,1], a[1,1])
        
        @printf("X / Y / przybliżenie / błąd względny / błąd bezwględny")
        
        minBW, maxBW = 10, 0
        minBB, maxBB = 10, 0
        
        ## Dodatkowy output potrzebny we sprawozdaniu
        
        for i in 1:len
            g = a[2,1]*matX[2, i] + a[1,1]
            
            bw = abs(valuesY[i] - g)
            bb = abs((valuesY[i] - g)/ valuesY[i])
            
            minBW = bw < minBW ? bw : minBW
            maxBW = bw > maxBW ? bw : maxBW
            
            minBB = bb < minBB ? bb : minBB
            maxBB = bb > maxBB ? bb : maxBB
            
            @printf("%d: %.10f %.10f %.10f %.10f %.10f\n",i, argumentsX[i], valuesY[i], g, bw, bb)
        end
        
        @printf("maxBW: %.10f, minBW: %.10f \nmaxBB: %.10f, minBB: %.10f\n", maxBW, minBW, maxBB, minBB)
        # wyniki sparsowane do wklejenia do LaTeXa
        @printf("%.10f & %.10f & %.3e & %.3e & %.3e & %.3e\n", a[2,1], a[1,1], minBB, maxBB, minBW, maxBW)
    end
    return a[2,1], a[1,1] 
end
#	====================================================================
#							Generatory zbiorów
#   ====================================================================
# Przedziały równoodległe, o punkcie początkowym begX, końcowym endX i kroku step
function inputgen(begX, endX, step)
    res = []
    x = begX
    while(x < endX)
        push!(res, x)
        x = x + step
    end
    
    return res
end

# zbiór X dla pierwszej metody (zdefiniowany we sprawozdaniu)
function definputs()
	res = []
	for i in 3:5
		for j in 1:9
			push!(res, j*10^i)
		end
	end
	return res
end
#	====================================================================
#							Regresja liniowa
#   ====================================================================
function compute(inputs; verbose = false, chartName = "", chartTitle = "")
    valuesY = []
    argumentsX = []
    chartY = []
    
    for i in 1:length(inputs)
        push!(chartY, gamma(inputs[i]) -eulerConstant)
        push!(valuesY, log2(chartY[i]))
        push!(argumentsX, log2(inputs[i]))
    end
    
    a, b = linear_regression(argumentsX, valuesY, verbose)
    
    c = 2^b
    d = -a
    
    
    
    if verbose       
        @printf("c = %.20f\nd = %.20f\n", c, d)   
    end
    
    # generator wykresu
    if chartName != ""
		f(x) = c*x^(-d)
    
		resPlot = plot(LinRange(inputs[1], inputs[length(inputs)], 100), f , color="red", linewidth=2.0, #lineslista1/Untitled3tyle=:dash, 
        title=chartTitle, legend=:outertopleft, xlabel="\$x\$", label ="\$ cn^{-d} \$")

        Plots.scatter!(inputs, chartY , color="blue", label="\$ \\gamma_n - \\gamma \$" )
        
        savefig(resPlot, chartName)
    end
    
    return c, d
end
#	====================================================================
#							Druga Metoda
#   ====================================================================
function error(n)
	return gamma(n) - eulerConstant
end

function convergenceExponent(n)
	return log2(error(n)/error(n-1))/log2(error(n-1)/error(n-2))
end

# druga metoda
function compute2(inputs; verbose = false, chartName = "", chartTitle = "")
    valuesY = []
    
    # inicjacja tablic z wartościami
    for i in 1:length(inputs)
        push!(valuesY, convergenceExponent(inputs[i]))
    end 
    
    
    # rysowanie wykresu
    if chartName != ""

        rplot = Plots.scatter(inputs, valuesY , color="blue", title=chartTitle, legend=:outertopleft, label ="\$ c \$")
        
        f(x) = 1
        
        plot!(LinRange(inputs[1], inputs[length(inputs)], 100), f , color="red", label="\$ f(x) = 1 \$" )
    
		savefig(rplot, chartName)
    end
end


#	====================================================================
#							Porównanie
#   ====================================================================
function compare(inputs; verbose = false, chartName = "", chartTitle = "")
    valuesY1 = []
    valuesY2 = []
    
    # inicjacja tablic z wartościami
    for i in 1:length(inputs)
        c, d = compute([inputs[i]-2, inputs[i]-1, inputs[i]])
        push!(valuesY1, d)
        push!(valuesY2, convergenceExponent(inputs[i]))
    end 
    
    
    # wykres
    if chartName != ""

        rplot = Plots.scatter(inputs, valuesY1,  color="blue", title=chartTitle, legend=:outertopleft, label ="Metoda 1")
        Plots.scatter!(inputs, valuesY2, color="yellow", label ="Metoda 2")
        
        f(x) = 1
        
        plot!(LinRange(inputs[1], inputs[length(inputs)], 100), f , color="red", label="\$ f(x) = 1 \$" )
        
        savefig(rplot, chartName)
    end
end

#	====================================================================
#							Wykresy i wyniki
#   ====================================================================

# wybor trybu dla programu
mode = 4

# tryb 1 - metoda 1
# tryb 2 - metoda 2
# tryb 3 - porównanie metod
# tryb 4 - wszystkie tryby jednocześnie

## Metoda 1
if mode == 1 || mode == 4
	compute(definputs(); 				 verbose = true, chartName = "wykresy/chart1_1.png", chartTitle="Metoda 1 - Zbiór \$ X \$")
	compute(inputgen(10^2, 10^2+100, 3); verbose = true)
	compute(inputgen(10^3, 10^3+100, 3); verbose = true)
	compute(inputgen(10^4, 10^4+100, 3); verbose = true)
	compute(inputgen(10^5, 10^5+100, 3); verbose = true)
	compute(inputgen(10^6, 10^6+100, 3); verbose = true, chartName = "wykresy/chart1_6.png", chartTitle="Metoda 1 - Zbiór \$ X_{6} \$")
end

## Metoda 2
if mode == 2 || mode == 4
	compute2(inputgen(10^6, 10^6+1000, 20); verbose = true, chartName = "wykresy/chart2_3.png", chartTitle="Metoda 2 - Zbiór \$ A_{3} \$")
	compute2(inputgen(10^4, 10^4+1000, 20); verbose = true, chartName = "wykresy/chart2_2.png", chartTitle="Metoda 2 - Zbiór \$ A_{2} \$")
	compute2(inputgen(10^3, 10^3+1000, 20); verbose = true, chartName = "wykresy/chart2_1.png", chartTitle="Metoda 2 - Zbiór \$ A_{1} \$")
end

## Porównanie
if mode == 3 || mode == 4
	compare(inputgen(10^6, 10^6+1000, 20); verbose = true, chartName = "wykresy/chart3_3.png", chartTitle="Porównanie Metod - Zbiór \$ A_{3} \$")
	compare(inputgen(10^4, 10^4+1000, 20); verbose = true, chartName = "wykresy/chart3_2.png", chartTitle="Porównanie Metod - Zbiór \$ A_{2} \$")
	compare(inputgen(10^3, 10^3+1000, 20); verbose = true, chartName = "wykresy/chart3_1.png", chartTitle="Porównanie Metod - Zbiór \$ A_{1} \$")
end



