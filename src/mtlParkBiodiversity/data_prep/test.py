import importlib.util

# For illustrative purposes.
package_name = 'mtlParkBiodiversity'

if importlib.util.find_spec(package_name) is None:
    print(package_name +" is not installed")

import mtlParkBiodiversity