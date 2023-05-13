import tkinter as tk
from PIL import ImageTk, Image
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np

# Constants
AU = 149.6e6  			 	# Astronomical Unit in kilometers
SolarConstant = 1361.0   	# Solar constant in W/m^2
sigma = 5.67e-8          	# Stefen Boltzmann Costant

# These are just for reference
CO2Concentration = 0.0004   # Concentration of CO2 in Earth's atmosphere (400 ppm)
CH4Concentration = 1.8e-6	# Concentration of CH4 in Earth's atmosphere (1.8 ppm)


# Calculate the surface temperature of the exoplanet
def calculate_surface_temperature(sigma, flux, co2, ch4, radius, pressure, density, heat):
	temp = ((flux * (1.0 - co2 - ch4)) / (4.0 * radius**2 * sigma)) ** 0.25
	temperature = temp * (1.0 + ((co2 + ch4) * pressure * 1000) / (density * heat))
	return temperature

# Handle button click event
def calculate_temperature():
    try:
        distance = float(distance_entry.get())
        radius = float(radius_entry.get())
        pressure = float(pressure_entry.get())
        density = float(density_entry.get())
        heat = float(specificHeat_entry.get())
        co2 = float(co2_entry.get())
        ch4 = float(ch4_entry.get())

        # Calculate the surface temperature
        F = SolarConstant / (distance **2)
        CO2Parameter = co2 / 1000000 * 5.35 # Scaling factor for CO2 greenhouse effect (Earth-like climate) - https://www.ipcc.ch/report/ar5/wg1/
        CH4Parameter = ch4 / 1000000 * 0.036 # Scaling factor for CH4 greenhouse effect (Earth-like climate) - https://www.ipcc.ch/report/ar5/wg1/
        temperature = calculate_surface_temperature(sigma, F, CO2Parameter, CH4Parameter, radius, pressure, density, heat)

        # Update the result label
        result_label.configure(text="Surface Temperature: {:.2f} °C".format(temperature-273.15))
        # Update the circle size and color
        circle_size = min(int(radius * 10), 95)
        circle_color = rgbtohex(get_rgb(temperature))
        canvas.coords(circle, 100 - circle_size, 100 - circle_size, 100 + circle_size, 100 + circle_size)
        canvas.itemconfigure(circle, width=circle_size, fill=circle_color) 
        
        plot_comparison()
              
    except ValueError:
        result_label.configure(text="Invalid input. Please enter numeric values.")

# Get color based on temperature
        
def get_rgb(temperature):
    # You can customize the color mapping based on your preference
	if temperature < 200:
		r = 0 
		g = 0
		b = 255
	elif temperature < 273:
        # Blue to cyan gradient
		r = 0
		g = int(255 * (temperature - 200) / 73)
		b = 255
	elif temperature < 400:
        # Cyan to green gradient
		r = 0
		g = 255
		b = int(255 * (400 - temperature) / 127)
	elif temperature < 600:
        # Green to yellow gradient
		r = int(255 * (temperature - 400) / 199)
		g = 255
		b = 0
	else:
        # Yellow to red gradient
		r = 255
		g = max(int(255 * (1000 - temperature) / 400), 0)
		b = 0

	return (r, g, b)

    
def rgbtohex(rgb):
	r, g, b, = rgb
	return f"#{r:02x}{g:02x}{b:02x}"
	
	
# Bar graph

def plot_comparison():
    parameters = ['D','R', 'P', 'Rho', 'Q', 'CO2', 'CH4']
    exoplanet_values = [float(distance_entry.get()),float(radius_entry.get()), float(pressure_entry.get()), float(density_entry.get()),
                        float(specificHeat_entry.get()), float(co2_entry.get()), float(ch4_entry.get())]
    earth_values = [1.0 ,1.0, 100.0, 1.2, 700.0, 400.0, 1.8]

    fig, ax = plt.subplots(figsize=(3,4))
    bar_width = 0.3  # Adjust the bar width

    # Set the positions of the bars
    x_pos = np.arange(len(parameters))

    ax.bar(x_pos, exoplanet_values, width=bar_width,label='Exoplanet')
    ax.bar(x_pos + bar_width, earth_values, width=bar_width, label='Earth')

    ax.set_ylabel('Values')
    ax.set_title('Exoplanet vs Earth')
    ax.set_xticks(x_pos + bar_width / 2)
    ax.set_xticklabels(parameters, rotation=45, ha='right')  # Rotate x-axis labels by 45 degrees
    ax.legend()

    # Clear previous content in the bar_canvas
    bar_canvas.delete('all')

    # Create a FigureCanvasTkAgg instance and display it in the bar_canvas
    canvas2 = FigureCanvasTkAgg(fig, master=bar_canvas)
    canvas2.draw()
    canvas2.get_tk_widget().grid(row=0, column=0)

    # Add a toolbar
    toolbar = NavigationToolbar2Tk(canvas2, bar_canvas)
    toolbar.update()
    canvas2.get_tk_widget().grid(row=0, column=0)
	
	
# Create the GUI
window = tk.Tk()
window.title("Exoplanet Temperature Simulator")

# Result label
result_label = tk.Label(window, text=" ", font=("Arial", 14))
result_label.grid(row=9, column=0, columnspan=2)



# Additional Notes Box
T = tk.Text(window, height=10, width=30)
T.grid(row=0, column=0, columnspan=2, sticky=tk.W+tk.E)
quote = """This simulator gives a rough estimate for an exoplanet 
surface temperature.
The host star is assumed to be sun-like (1361 W/m^2 
irradiance).
The chemical composition is only used to dampen the flux 
and increase the weight of pressure.
Earth parameters in given units for reference: P = 100,
Density = 1.2, Sp.Heat=700, CO2 = 400, CH4 = 1.8

Marco Leonardi - 05/2023 - Unibo"""
T.insert(tk.END, quote)
T.config(state="disabled")


# Distance input
distance_label = tk.Label(window, text="Distance (AU):")
distance_label.grid(row=1, column=0)
distance_entry = tk.Entry(window)
distance_entry.grid(row=1, column=1)

# Radius input
radius_label = tk.Label(window, text="Radius (Earth radii):")
radius_label.grid(row=2, column=0)
radius_entry = tk.Entry(window)
radius_entry.grid(row=2, column=1)

# Pressure input
pressure_label = tk.Label(window, text="Pressure (kPa)")
pressure_label.grid(row=3, column=0)
pressure_entry = tk.Entry(window)
pressure_entry.grid(row=3, column=1)

# Density input
density_label = tk.Label(window, text="Density (kg/m^-3)")
density_label.grid(row=4, column=0)
density_entry = tk.Entry(window)
density_entry.grid(row=4, column=1)

# Heat input
heat_label = tk.Label(window, text="Specific Heat (J/Kg*K)")
heat_label.grid(row=5, column=0)
specificHeat_entry = tk.Entry(window)
specificHeat_entry.grid(row=5, column=1)

# CO2 input
co2_label = tk.Label(window, text="Co2 Concentration (ppm \\ mg/L)")
co2_label.grid(row=6, column=0)
co2_entry = tk.Entry(window)
co2_entry.grid(row=6, column=1)

# CH4 input
ch4_label = tk.Label(window, text="Methane Concentration (ppm \\ mg/L)")
ch4_label.grid(row=7, column=0)
ch4_entry = tk.Entry(window)
ch4_entry.grid(row=7, column=1)

# Calculate button
calculate_button = tk.Button(window, text="Calculate", command=calculate_temperature)
calculate_button.grid(row=8, column=0, columnspan=2)

# Plot button
#plot_button = tk.Button(window, text="Plot", command=plot_comparison)
#plot_button.grid(row=8, column=3)

# Circle widget
canvas = tk.Canvas(window, width=250, height=250)
circle = canvas.create_oval(100, 100, 100, 100, outline="")
canvas.grid(row=15, column=0, columnspan=3)

# Bar graph canvas
bar_canvas = tk.Canvas(window, width=250, height=250)
bar_canvas.grid(row=15, column=3, columnspan = 3)


# Start the GUI event loop
window.mainloop()
