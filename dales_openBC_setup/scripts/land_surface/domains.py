
# x0, y0 are in RD coordinates

eindhoven = {'expname': 'eindhoven',
             'x0'     : 152000,
             'y0'     : 375000,
             'itot'   : 440,
             'jtot'   : 320,
             'dx'     : 50,
             'dy'     : 50,
             'nprocx' : 4,
             'nprocy' : 4
            }

test       = {'expname': 'test',
             'x0'     : 912500,
             'y0'     : 940000,
             'itot'   : 64,
             'jtot'   : 64,
             'dx'     : 400,
             'dy'     : 400,
             'nprocx' : 2,
             'nprocy' : 2
            }

ruisdael   = {'expname': 'ruisdael',
             'x0'     : 910000,
             'y0'     : 940000,
             'itot'   : 864,
             'jtot'   : 576,
             'dx'     : 200,
             'dy'     : 200,
             'nprocx' : 12,
             'nprocy' : 12
            }

smalldomain = {'expname': 'smalldomain',
            'x0' : 155000,
            'y0' : 386000,
            'itot' : 128,
            'jtot' : 128,
            'dx' : 50,
            'dy' : 50,
            'nprocx' : 4,
            'nprocy' : 4
            }

veluwe = {'expname': 'veluwe',
            'x0' : 174000,
            'y0' : 448500,
            'itot' : 400,
            'jtot' : 200,
            'dx' : 50,
            'dy' : 50,
            'nprocx' : 8,
            'nprocy' : 4
            }

veluwe_small = {'expname': 'veluwe_small',
            'x0' : 186000,
            'y0' : 448500,
            'itot' : 160,
            'jtot' : 80,
            'dx' : 50,
            'dy' : 50,
            'nprocx' : 4,
            'nprocy' : 2
            }

gouda = {'expname': 'gouda',
            'x0' : 100942,
            'y0' : 446028,
            'itot' : 128,
            'jtot' : 128,
            'dx' : 625,
            'dy' : 625,
            'nprocx' : 2,
            'nprocy' : 2
            }

ruisdael = {'expname': 'ruisdael', # large Ruisdael domain
            'x0' : -98927.04501275037,
            'y0' : 296195.8715821294,
            'itot' : 768,
            'jtot' : 512,
            'dx' : 625,
            'dy' : 625,
            'nprocx' : 12,
            'nprocy' : 16,
            }
            
rotterdam_nest = {'expname': 'rotterdam_nest', # rotterdam_nested domain
            'x0' : 79960.0, 
            'y0' : 426325.0,
            'itot' : 512,
            'jtot' : 512,
            'dx' : 50,
            'dy' : 50,
            'nprocx' : 16,
            'nprocy' : 16,
            }
            
            
            
rotterdam_new_coarse = {'expname': 'rotterdam_new_coarse', # rotterdam_new_coarse domain
            'x0' : 76125.0, 
            'y0' : 396285.0,
            'itot' : 256,
            'jtot' : 256,
            'dx' : 400,
            'dy' : 400,
            'nprocx' : 24,
            'nprocy' : 16,
            }
            
co2_ruisdael_new_fine = {'expname': 'co2_ruisdael_new_fine', # rotterdam_new_coarse domain
            'x0' : 52440.0, 
            'y0' : 419115.0,
            'itot' : 1024,
            'jtot' : 512,
            'dx' : 100,
            'dy' : 100,
            'nprocx' : 32,
            'nprocy' : 16,
            }
            
            
domains = {'eindhoven': eindhoven,
           'eindhoven_small': smalldomain,
           'test'     : test,
           'ruisdael' : ruisdael,
           'veluwe' : veluwe,
           'veluwe_small' : veluwe_small,
           'gouda' : gouda,
           'ruisdael' : ruisdael,
           'rotterdam_nest' : rotterdam_nest,
           'rotterdam_new_coarse' : rotterdam_new_coarse,
           'co2_ruisdael_new_fine' : co2_ruisdael_new_fine
           }
