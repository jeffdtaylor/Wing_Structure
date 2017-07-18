import optix
import jeff

model = optix.objective_model(jeff.script, max_processes = 8)
#settings = optix.settings.load('settings.json')
settings = optix.settings()
settings.add_variable('weight', 3000.0)
settings.add_variable('wingspan', 89.0)

optix.optimize(model, settings)

