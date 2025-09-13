from .park_metrics import park_metrics

def process_all_metrics(force = False, test = False, limit = None):

    # Run park metrics
    park_metrics(force = force, test = test, limit = limit)

    #
