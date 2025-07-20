import tkinter as tk
from tkinter import ttk, messagebox
import math

class SolarCollectorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Solar Flat Plate Collector Efficiency Calculator")
        self.root.geometry("800x700") # Adjust window size as needed

        self.create_widgets()

    def create_widgets(self):
        # Create a notebook (tabbed interface) for better organization
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(pady=10, expand=True, fill="both")

        # --- Input Frame ---
        self.input_frame = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(self.input_frame, text="Input Parameters")
        self.create_input_widgets(self.input_frame)

        # --- Results Frame ---
        self.results_frame = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(self.results_frame, text="Calculation Results")
        self.create_results_widgets(self.results_frame)

    def create_input_widgets(self, parent_frame):
        # Helper function to create a labeled input field
        def create_input_row(frame, label_text, default_value=""):
            row = ttk.Frame(frame)
            row.pack(pady=5, fill="x")
            label = ttk.Label(row, text=label_text, width=30)
            label.pack(side="left", padx=5)
            entry = ttk.Entry(row, width=40)
            entry.insert(0, str(default_value))
            entry.pack(side="right", expand=True, fill="x")
            return entry

        ttk.Label(parent_frame, text="Enter Collector and Location Details:", font=("Arial", 12, "bold")).pack(pady=10)

        self.beta_entry = create_input_row(parent_frame, "Tilt Angle (degrees):", 30) # Example default
        self.n_entry = create_input_row(parent_frame, "Day of the Year (1-365):", 172) # Example default (June 21)
        self.phi_entry = create_input_row(parent_frame, "Latitude of Location (degrees):", 27.6) # Example default (approx Kota)
        self.longitude_entry = create_input_row(parent_frame, "Longitude of Location (degrees):", 75.8) # Example default (approx Kota)

        ttk.Label(parent_frame, text="Time of Day (24-hour format):", font=("Arial", 10, "bold")).pack(pady=5)
        self.hours_entry = create_input_row(parent_frame, "Hours:", 12)
        self.minutes_entry = create_input_row(parent_frame, "Minutes:", 30)

        ttk.Label(parent_frame, text="Flat Plate Collector Dimensions:", font=("Arial", 10, "bold")).pack(pady=5)
        self.length_entry = create_input_row(parent_frame, "Length (m):", 2.0)
        self.breadth_entry = create_input_row(parent_frame, "Breadth (m):", 1.0)
        self.height_entry = create_input_row(parent_frame, "Height (m) (for Us calculation, typically insulation thickness):", 0.05) # Assuming insulation thickness for H in Us
        self.num_glass_covers_entry = create_input_row(parent_frame, "Number of Glass Covers (M):", 1)

        self.calculate_button = ttk.Button(parent_frame, text="Calculate Efficiency", command=self.perform_calculation)
        self.calculate_button.pack(pady=20)

    def create_results_widgets(self, parent_frame):
        ttk.Label(parent_frame, text="Calculation Results:", font=("Arial", 12, "bold")).pack(pady=10)

        self.declination_label = ttk.Label(parent_frame, text="Declination Angle: ")
        self.declination_label.pack(anchor="w", pady=2)

        self.hour_angle_label = ttk.Label(parent_frame, text="Hour Angle: ")
        self.hour_angle_label.pack(anchor="w", pady=2)

        self.total_radiation_label = ttk.Label(parent_frame, text="Total Radiation (It): ")
        self.total_radiation_label.pack(anchor="w", pady=2)

        self.useful_heat_label = ttk.Label(parent_frame, text="Useful Heat Gain: ")
        self.useful_heat_label.pack(anchor="w", pady=2)

        self.efficiency_label = ttk.Label(parent_frame, text="Collector Efficiency: ")
        self.efficiency_label.pack(anchor="w", pady=2)

        # You can add more labels here for other intermediate results if desired
        self.ibn_label = ttk.Label(parent_frame, text="Normal Beam Radiation (Ibn): ")
        self.ibn_label.pack(anchor="w", pady=2)
        self.ib_label = ttk.Label(parent_frame, text="Normal Radiation (Ib): ")
        self.ib_label.pack(anchor="w", pady=2)
        self.id_label = ttk.Label(parent_frame, text="Diffused Radiation (Id): ")
        self.id_label.pack(anchor="w", pady=2)
        self.ig_label = ttk.Label(parent_frame, text="Global Radiation (Ig): ")
        self.ig_label.pack(anchor="w", pady=2)
        self.ul_label = ttk.Label(parent_frame, text="Overall Loss Coefficient (Ul): ")
        self.ul_label.pack(anchor="w", pady=2)
        self.s_label = ttk.Label(parent_frame, text="Absorbed Incident Radiation (S): ")
        self.s_label.pack(anchor="w", pady=2)


    def perform_calculation(self):
        try:
            # Get values from input fields
            beta = float(self.beta_entry.get())
            n = int(self.n_entry.get())
            phi = float(self.phi_entry.get())
            longitude = float(self.longitude_entry.get())
            hours = int(self.hours_entry.get())
            minutes = int(self.minutes_entry.get())
            L = float(self.length_entry.get())
            B = float(self.breadth_entry.get())
            H = float(self.height_entry.get()) # Used for Us, assumed insulation thickness
            M = int(self.num_glass_covers_entry.get())

            # --- Your Calculation Logic (from your provided code) ---
            # All angle conversions to radians [cite: 2]
            beta_r = math.radians(beta)
            phi_r = math.radians(phi)

            # Calculating declination angle delta [cite: 3]
            delta = (math.sin((284 + n) * 360 / 365)) * 23.45
            delta_r = math.radians(delta)

            # Calculating solar time or LAT from the time given [cite: 4]
            total_time = (hours * 60) + minutes
            B_val = (n - 81) * 360 / 365 # Renamed 'B' to 'B_val' to avoid conflict with collector breadth
            Br = math.radians(B_val)
            EOT = 9.87 * math.sin(Br * 2) - 7.53 * math.cos(Br) - 1.5 * math.sin(Br)
            x = (4 * (82.5 - longitude))
            LAT = total_time - x - EOT
            omega = (720 - LAT) * 15 / 60
            omega_r = math.radians(omega)

            # Cosine of zenith angle for beam radiation calculations
            cosazm = (math.sin(phi_r) * math.sin(delta_r)) + (math.cos(phi_r) * math.cos(delta_r) * math.cos(omega_r))

            # Check for valid cosazm to avoid math domain errors (e.g., during night)
            if cosazm <= 0:
                messagebox.showinfo("Result", "No direct solar radiation (Sun below horizon or at very low angle). Cannot calculate efficiency based on current model.")
                # Clear or set results to 0/N/A
                self.declination_label.config(text=f"Declination Angle: {delta:.2f} degrees")
                self.hour_angle_label.config(text=f"Hour Angle: {omega:.2f} degrees")
                self.total_radiation_label.config(text="Total Radiation (It): N/A")
                self.useful_heat_label.config(text="Useful Heat Gain: N/A")
                self.efficiency_label.config(text="Collector Efficiency: N/A")
                self.ibn_label.config(text="Normal Beam Radiation (Ibn): N/A")
                self.ib_label.config(text="Normal Radiation (Ib): N/A")
                self.id_label.config(text="Diffused Radiation (Id): N/A")
                self.ig_label.config(text="Global Radiation (Ig): N/A")
                self.ul_label.config(text="Overall Loss Coefficient (Ul): N/A")
                self.s_label.config(text="Absorbed Incident Radiation (S): N/A")
                return

            A = 1202
            B_rad_const = 0.141 # Renamed to avoid conflict with collector breadth
            C_rad_const = 0.103 # Renamed to avoid conflict with Tkinter 'c'
            Ibn = A * math.exp((-1) * B_rad_const / cosazm) # Normal beam radiation [cite: 4]
            Ib = Ibn * cosazm # Normal radiation [cite: 5]
            Id = C_rad_const * cosazm # Diffused radiation [cite: 5]
            Ig = Ib + Id # Global radiation [cite: 6]

            phiminusbeta = phi - beta
            phiminusbetar = math.radians(phiminusbeta)

            rb = (math.sin(delta_r) * math.sin(phiminusbetar) + math.cos(phiminusbetar) * math.cos(omega_r) * math.cos(delta_r)) / \
                 (math.sin(delta_r) * math.sin(phi_r) + (math.cos(phi_r) * math.cos(delta_r) * math.cos(omega_r)))
            rd = (1 + math.cos(beta_r)) / 2
            rho = 0.2
            rr = rho * (1 - math.cos(beta_r)) / 2
            It = (Ib * rb) + (Id * rd) + (Ig * rr) # Total radiation on a tilted surface [cite: 6]

            # Area of collector - assuming L*B is the aperture area
            area = L * B

            # Assuming standard values of emissivity and transmissivity
            emissivity = 0.95
            transmissivity = 0.90

            # Incident radiation that is absorbed by plate
            # The original code had a syntax error (Ib*rb(emissivity*transmissivity)). Corrected assuming multiplication.
            S = (Ib * rb * (emissivity * transmissivity)) + ((Id * rd) + (Ig * rr)) * (emissivity * transmissivity)

            # Assuming hw, thermal conductivity of insulation (ki)
            hw = 7.889
            ki = 0.04

            # Assuming 'delta' in Us and Ub is the insulation thickness. If H is insulation thickness, then use H.
            # Your original code used 'delta' (declination angle) for thickness, which is incorrect.
            # I am assuming H is the insulation thickness for the Ub and Us calculations.
            # If 'delta' in your original code for Us and Ub actually refers to an insulation thickness, please clarify.
            # For now, I'm using H for insulation thickness.
            # Also, Us calculation (L+B)*ki*H/(L*B*delta) seems to imply delta is insulation thickness for side/back loss.
            # Let's assume 'delta_insulation' is implied for Ub and Us, and H is that value.

            # Re-evaluating Us and Ub based on standard practices:
            # Us (Side Loss): For typical flat plate collectors, side losses are often approximated as a perimeter loss.
            # Ub (Back Loss): ki / insulation_thickness
            # Given your code:
            # Us = ((L+B)*ki*H)/(L*B*delta) -> If delta here means insulation thickness, then replace it with H.
            # Ub = ki/delta -> If delta here means insulation thickness, then replace it with H.
            # Let's use H as the insulation thickness for these calculations.
            # If the variable 'delta' in your Us and Ub calculation was *intended* to be declination angle, it's physically incorrect for heat transfer.
            # Assuming 'delta' for insulation thickness was a placeholder in your original code:
            
            if H <= 0:
                raise ValueError("Height (insulation thickness) must be greater than 0 for loss calculations.")

            Us = ((L + B) * ki * H) / (L * B * H) # If H is insulation thickness, and delta was meant to be insulation thickness.
            Us = (L+B)*ki / (L*B) # A simplified interpretation if H cancels out. Let's use your formula directly and assume H is thickness
            # Let's re-interpret the original Us and Ub more closely, assuming a typo where 'delta' was meant to be H (insulation thickness)
            # based on its context in denominator for conductivity.
            
            # Corrected Us and Ub assuming 'delta' in your original code for Us and Ub was meant to be 'H' (insulation thickness)
            Us_corrected = ((L + B) * ki * H) / (L * B * H) if H > 0 else 0 # Sides area / Collector Area * ki/thickness
            if H > 0:
                Ub_corrected = ki / H
            else:
                Ub_corrected = 0 # Or handle error appropriately

            # Assuming typical operating temperatures
            Tp = 373  # Collector Plate Temperature in Kelvin (100 C)
            Ta = 298  # Ambient Air Temperature in Kelvin (25 C)
            sigma = 5.67e-8  # Stefan-Boltzmann constant (W/m^2.K^4)  - Original code had 567000000 which is incorrect magnitude

            f = (1 - 0.04 * hw + 0.005 * hw * hw) * (1 + 0.091 * M)
            c = 365.9 * (1 - 0.00883 * beta + 0.0001298 * beta * beta) # Beta in degrees for this formula

            var1 = (Tp - Ta) / (M + f)

            # Ut calculation [cite: 8]
            # Be careful with the power operator 'pow' and the entire structure.
            # The original Ut formula looks like it has a potential formatting issue or is a very specific empirical correlation.
            # Let's break it down as it was given, but be aware empirical formulas can be very sensitive.
            # (sigma*(Tp*Tp - Ta*Ta)*(Tp - Ta)) should be sigma*(Tp^4 - Ta^4) for radiation.
            # Your formula is (sigma*(Tp*Tp - Ta*Ta)*(Tp - Ta)). This is not sigma*(Tp^4 - Ta^4).
            # The provided formula for Ut is complex and might be a specific empirical correlation.
            # I'll directly translate your provided formula for Ut.
            
            term1_denominator_part1 = M / (c / Tp) * (pow(var1, 0.33))
            term1_denominator = term1_denominator_part1 + (1 / hw)
            term1 = pow(term1_denominator, -1)

            # Radiation term (original code): (sigma*(Tp*Tp - Ta*Ta)*(Tp - Ta))/(1/(emissivity+0.005*M*(1-emissivity)))+((2*M+f-1)/emissivity)-M)
            # This looks like (sigma * (T_p^2 - T_a^2) * (T_p - T_a)) / Denominator
            # If it's meant to be Stefan-Boltzmann, it should be sigma * (T_p^4 - T_a^4)
            # I will use your exact formula from the source, but note this deviation from standard radiation terms.
            
            radiation_term_numerator = (sigma * (Tp*Tp - Ta*Ta) * (Tp - Ta))
            radiation_term_denominator_part1 = 1 / (emissivity + 0.005 * M * (1 - emissivity))
            radiation_term_denominator_part2 = ((2 * M + f - 1) / emissivity) - M
            radiation_term_denominator = radiation_term_denominator_part1 + radiation_term_denominator_part2
            
            # Avoid division by zero
            if radiation_term_denominator == 0:
                rad_term = float('inf') # Or handle as an error
            else:
                rad_term = radiation_term_numerator / radiation_term_denominator

            Ut = term1 + rad_term # [cite: 8]

            Ul = Ub_corrected + Us_corrected + Ut # Overall loss coefficient

            usefulheat = (S - Ul * (Tp - Ta))
            
            if It == 0:
                efficiency = 0.0 # Avoid division by zero if no radiation
            else:
                efficiency = usefulheat / It

            # Update results in GUI
            self.declination_label.config(text=f"Declination Angle: {delta:.2f} degrees")
            self.hour_angle_label.config(text=f"Hour Angle: {omega:.2f} degrees")
            self.total_radiation_label.config(text=f"Total Radiation (It): {It:.2f} W/sq.m")
            self.useful_heat_label.config(text=f"Useful Heat Gain: {usefulheat:.2f} W/sq.m")
            self.efficiency_label.config(text=f"Collector Efficiency: {efficiency*100:.2f}%")
            self.ibn_label.config(text=f"Normal Beam Radiation (Ibn): {Ibn:.2f} W/sq.m")
            self.ib_label.config(text=f"Normal Radiation (Ib): {Ib:.2f} W/sq.m")
            self.id_label.config(text=f"Diffused Radiation (Id): {Id:.2f} W/sq.m")
            self.ig_label.config(text=f"Global Radiation (Ig): {Ig:.2f} W/sq.m")
            self.ul_label.config(text=f"Overall Loss Coefficient (Ul): {Ul:.2f} W/m^2.K")
            self.s_label.config(text=f"Absorbed Incident Radiation (S): {S:.2f} W/sq.m")


            self.notebook.select(self.results_frame) # Switch to results tab

        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid input. Please ensure all fields are numeric and correct. Error: {e}")
        except Exception as e:
            messagebox.showerror("Calculation Error", f"An error occurred during calculation: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = SolarCollectorApp(root)
    root.mainloop()