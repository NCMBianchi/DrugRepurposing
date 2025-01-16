"""
FORMS MODULE: DEFINES THE USER INPUTS, IN EACH HTML FORM
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import multiprocessing
from flask_wtf import FlaskForm
from flask import session
from wtforms import StringField, SubmitField, SelectField, BooleanField, IntegerField, DecimalField
from wtforms.validators import DataRequired, Length, NumberRange, Optional

class user_input(FlaskForm):

    # date toggle (1 if flagged, or 0 if not)
    date_t = BooleanField('Date Override', default=False)

    # date override (Day, Month, Year)
    date_OR_day = IntegerField('Override Day', validators=[NumberRange(min=1, max=31)], default=31)
    date_OR_month = IntegerField('Override Month', validators=[NumberRange(min=1, max=12)], default=7)
    date_OR_year = IntegerField('Override Year', validators=[NumberRange(min=2020, max=2050)], default=2024)
    # eventually change the upper bound for the year

    # disease URI/ID
    disease_URI = StringField('Insert the \'MONDO:\' URI for the disease of interest\n(e.g. MONDO:0007739 for Huntington Disease)', 
                               validators=[DataRequired(), Length(min=5, max=20)], default='MONDO:0007739')
    
    # degree of distance (drop down) [denugger mode adds 1 as an option: TO REMOVE]
    #deg_of_dist = SelectField('Degrees of Distance', choices=[(2, '2'), (3, '3')], default=3, coerce=int)
    deg_of_dist = SelectField('Degrees of Distance', choices=[(1, '1'), (2, '2'), (3, '3')], default=3, coerce=int)

    # minimum similarity threshold
    inp_minimum_sim = DecimalField('Minimum Drug Similarity', validators=[NumberRange(min=0.0, max=1.0)], places=2, default=0.5, render_kw={"id": "inp_minimum_sim"})

    # negative samples toggle (1 if flagged, or 0 if not)
    ns_toggle = BooleanField('Negative Samples', default=True)

    # drug similarity threshold
    sim_t = DecimalField('Drug Similarity Threshold', validators=[NumberRange(min=0.5, max=1.0)], places=2, default=0.9, render_kw={"id": "sim_t"})
    # fix to avoid only 0 and 1

    # number of cores (from 1 to num_cores)
    num_cores = multiprocessing.cpu_count()
    n_cores = IntegerField('Number of CPU Cores', validators=[NumberRange(min=1, max=num_cores)], default=1)

    # mode of operation (drop down)
    po_mode = SelectField('Mode of Operation', choices=[('ultralight', 'Ultralight'), 
                                                        ('light', 'Light'), 
                                                        ('full', 'Full')],
                          default='light')

    # machine learning seed
    ML_seed = StringField('ML Seed', validators=[Optional()], default='random')

    # submit button
    submit = SubmitField('LAUNCH RUN')

    def validate_ML_seed(self, field):
        """
        Custom validation for ML Seed to accept integers or the word 'random'.
        """
        if field.data and field.data.lower() == 'random':
            field.data = 'random'
        elif field.data:
            try:
                field.data = int(field.data)
            except ValueError:
                raise ValueError("ML Seed must be an integer or 'Random'")