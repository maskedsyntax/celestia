import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def visualize_simulation(csv_file):
    # Load the data
    # Format: time, id, mass, x, y, z, vx, vy, vz
    cols = ['time', 'id', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']
    try:
        data = pd.read_csv(csv_file, names=cols)
    except FileNotFoundError:
        print(f"Error: {csv_file} not found. Run the simulation first.")
        return

    # Set up the plot
    fig = plt.figure(figsize=(10, 10), facecolor='black')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    ax.grid(False)
    ax.w_xaxis.pane.fill = False
    ax.w_yaxis.pane.fill = False
    ax.w_zaxis.pane.fill = False
    
    # Get unique time steps
    times = data['time'].unique()
    
    def update(frame):
        ax.clear()
        # Get data for current time
        current_data = data[data['time'] == times[frame]]
        
        # Plot particles
        # Scale size by mass (log scale for visibility)
        sizes = np.log10(current_data['mass'] + 1) * 10
        
        ax.scatter(current_data['x'], current_data['y'], current_data['z'], 
                   s=sizes, c='white', alpha=0.6, edgecolors='none')
        
        # Consistent bounds
        limit = 150
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_zlim(-limit, limit)
        
        ax.set_title(f'N-Body Simulation: Time {times[frame]:.3f}', color='white')
        ax.axis('off')
        return fig,

    print(f"Generating animation from {len(times)} frames...")
    ani = animation.FuncAnimation(fig, update, frames=len(times), interval=50)
    
    # Save the animation
    output_file = 'galaxy_simulation.gif'
    print(f"Saving to {output_file}...")
    ani.save(output_file, writer='pillow')
    print("Done!")

if __name__ == "__main__":
    visualize_simulation('galaxy_step.csv')
