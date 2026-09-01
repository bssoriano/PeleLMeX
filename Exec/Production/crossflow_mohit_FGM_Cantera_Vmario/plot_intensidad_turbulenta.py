#!/usr/bin/env python3
import yt
import numpy as np
import matplotlib.pyplot as plt

def main():
    # 1. Cargar el plotfile
    plt_file = "plt133981"  # Cambia esto por tu archivo objetivo
    try:
        ds = yt.load(plt_file)
    except Exception as e:
        print(f"Error al cargar el archivo: {e}")
        return

    # 2. Definir los parámetros físicos del chorro
    U_c = 2.0  # Velocidad característica del chorro (v_jet = 3.0 m/s)

    # 3. Crear una función para el campo derivado de Intensidad Turbulenta
    # Nota: Como es un plotfile individual, calcularemos la intensidad respecto
    # a la magnitud de la fluctuación espacial local o usando la velocidad local.
    # Aquí definimos la componente de intensidad turbulenta basándonos en la velocidad en Y:
    
    def _turbulent_intensity_y(field, data):
        # Obtenemos la velocidad en Y
        v_y = data[("boxlib", "y_velocity")].v
        
        # En una simulación instantánea, una aproximación común para u_rms local
        # sin promedios temporales previos es la desviación respecto al promedio del plano,
        # o simplemente evaluar la energía cinética de los componentes si posees u_rms guardado.
        # Si PeleLMeX no guardó 'u_rms' directamente, calculamos la fluctuación respecto a U_c:
        v_mean = np.mean(v_y) 
        u_rms = np.abs(v_y - v_mean) # Fluctuación local instantánea
        
        return u_rms / U_c

    # Registrar el nuevo campo en yt
    ds.add_field(
        name=("boxlib", "intensidad_turbulenta"),
        function=_turbulent_intensity_y,
        sampling_type="cell",
        units=""
    )

    print("Campo de 'intensidad_turbulenta' registrado exitosamente.")

    # 4. Crear un gráfico de contornos (SlicePlot) en el plano Z medio
    # Asumiendo que es un dominio 2D o queremos cortar en el centro de Z en 3D
    z_center = (ds.domain_left_edge[2] + ds.domain_right_edge[2]) / 2.0 if ds.dimensionality == 3 else 0.0
    
    # Generar el rebanado (Slice) perpendicular al eje Z para ver el perfil XY (o el plano del jet)
    # Si tu simulación es en un eje distinto, cambia 'z' por el eje normal al plano que deseas ver.
    slc = yt.SlicePlot(ds, 'z', ("boxlib", "intensidad_turbulenta"))
    
    # Ajustar parámetros estéticos del contorno
    slc.set_cmap(("boxlib", "intensidad_turbulenta"), "viridis") # Mapa de color suave para turbulencia
    slc.set_log(("boxlib", "intensidad_turbulenta"), False)     # Escala lineal para intensidades (%)
    slc.set_colorbar_label(("boxlib", "intensidad_turbulenta"), "Intensidad Turbulenta $I = u'_{rms} / U_c$")
    
    # Dibujar líneas de contorno explícitas sobre el mapa de colores
    # Esto añade las líneas de nivel (isolíneas)
    slc.annotate_contour(("boxlib", "intensidad_turbulenta"), levels=5, clim=(0.01, 0.20), plot_args={"colors": "white", "linewidths": 0.5})

    # Guardar la imagen del contorno
    output_image = "contorno_intensidad_turbulenta.png"
    slc.save(output_image)
    print(f"Gráfico de contornos guardado como: {output_image}")

if __name__ == "__main__":
    main()