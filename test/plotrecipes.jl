using GeoData, Test, Dates, Plots

ga2 = GeoArray(ones(91) * (-25:15)', (X(0.0:4.0:360.0), Y(-25.0:1.0:15.0), ); name=:Test)
ga3 = GeoArray(rand(10, 41, 91), (Z(100:100:1000), Y(-20.0:1.0:20.0), X(0.0:4.0:360.0)))
ga4ti = GeoArray(
    rand(10, 41, 91, 4), 
    (Z(100:100:1000), Y(-20.0:1.0:20.0), X(0.0:4.0:360.0), Ti(Date(2001):Year(1):Date(2004)))
)
ga4x = GeoArray(
    rand(10, 41, 91, 4), 
    (Z(100:100:1000), Y(-20.0:1.0:20.0), X(0.0:4.0:360.0), X())
)

plot(ga2)
plot(ga3[Y(At(0.0))])
plot(ga3[X(At(180.0))])
# Line plot with Z on the vertical axis
plot(ga3[X(At(0.0)), Y(At(0.0))])
# DD fallback line plot with Z as key (not great really)
plot(ga4x[X(At(0.0)), Y(At(0.0))])
# DD fallback heatmap with Z on Y axis
heatmap(ga4x[X(At(0.0)), Y(At(0.0))])
# Cant plot 4d
@test_throws ErrorException plot(ga4x)
# 3d plot by NoIndex X dim
plot(ga4x[Y(1)])
# 3d plot by Ti dim
plot(ga4ti[Z(1)])

# DD fallback
contour(ga2)
