import matplotlib

SEP_LITERAL = "="*50 + "\n"
PLT_SHAPES = ['o', 'p', "X", "s", 'p', 'P',
              '*', 'v', '^', '<', '>', '+', "x", "h"]

colormap = lambda name='Accent': matplotlib.cm.get_cmap(name)


def uv_sorted(e): return (e[0], e[1]) if e[0] < e[1] else (e[1], e[0])
