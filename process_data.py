def check_unit_system_consistency(dataframe):
    return True

def check_data_type(data,type):
    if type(data) == type:
        return True
    else:
        return False

def check_data_value(data,value):
    if data == value:
        return True
    else:
        return False

def define_plot_colors():
    from itertools import cycle
    colors = cycle(
        [
            'k',
            'b',
            'r',
            'g',
            'y',
            'c',
            'm'
        ]
    )
    return colors

def define_plot_lines():
    from itertools import cycle
    lines = cycle(
        [
            '-',
            '--',
            '-.',
            ':'
        ]
    )
    return lines

def define_plot_markers():
    from itertools import cycle
    markers = cycle(
        [
            '.',
            'o',
            'v',
            '^',
            '<',
            '>',
            '1',
            '2',
            '3',
            '4',
            '8',
            's',
            'p',
            'P',
            '*',
            'h',
            'H',
            '+',
            'x',
            'X',
            'D',
            'd'
        ]
    )
    return markers
