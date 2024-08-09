"""
FORMS MODULE: DEFINES THE USER INPUTS, IN EACH HTML FORM
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

from flask_wtf import FlaskForm
from flask import session
from wtforms import StringField, SubmitField, SelectField
from wtforms.validators import DataRequired, Length, Email, EqualTo


class user_input(FlaskForm):

    date_t = # 1 if flagged, or 0 if not (default = 0)

    date_OR_day = # optional, or necessary if month or year are filled 
    date_OR_month = # optional, or necessary if day or year are filled
    date_OR_year = # optional, or necessary if day or month are filled
    if date_OR_day and date_OR_month and date_OR_year:
        date_OR = datetime.date(date_OR_year, date_OR_month, date_OR_day)
    else:
        date_OR = datetime.date(2024, 7, 31)

    disease_URI = StringField('Insert the disease \'MONDO:\' URI (Monarch Initiative\'s own identifier) of the disease of interest (e.g. MONDO:0007739 for Huntington Disease) since OMIMs no longer are supported as query seeds.', validators =[DataRequired(), Length(min=5, max=20)])
    
    deg_of_dist = # 2 or 3, from a drop menu (default = 3)

    inp_minimum_sim = #

    ns_toggle = # 1 if flagged, or 0 if not (default = 1)

    sim_t = # (default = 0.9)

    n_cores = # 1 to num_cores
    
    po_mode = # from menu: 'ultralight', 'light', 'full'
    
    ML_seed = # integer as input from 1 to 1000000, (default = 'random')

    submit = SubmitField('LAUNCH RUN') 