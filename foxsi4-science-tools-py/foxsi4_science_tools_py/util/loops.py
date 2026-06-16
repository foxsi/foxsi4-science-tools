""" Functions to help with loops of all kinds. """

def progress_printer(counting, on, of_total, every_event=None):
    """Helper function to print out the progress of a loop. 
    
    Parameters
    ----------
    counting : `str`
        A helpful string describing what is being counted. 
    
    on : `int`
        The iteration the progress is on. The final value should be 
        expected to be equal to the `of_total` value for a new line 
        character to be added to the final statement (if `every_event` 
        is `None`).
    
    of_total : `int`
        The total number of iterations.
    
    every_event : `None` or `int`
        Instead of printing out every iteration, this number sets the 
        increment of each iteration. E.g., if `every_event=1_000` then 
        the progress will be printed out every 1_000 iterations. This 
        means the final printed percent may not be '100%'.
    """
    if every_event is None:
        end_line = "" if on<of_total else "\n"
        print(f"\r{counting} {on}/{of_total} ({((on)/of_total)*100:.1f}%)", end=end_line)
    else:
        if on%every_event==0:
            _total_loops = int(of_total/every_event) # find the total number of whole prints we need
            _on_loop = int(on/every_event) # find the current loop we're on
            end_line = "" if _on_loop<_total_loops else "\n"
            print(f"\r{counting} (print inc. {every_event}) {on}/{of_total} ({((on)/of_total)*100:.1f}%)", end=end_line)
