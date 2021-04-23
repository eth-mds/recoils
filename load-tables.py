import data_warehouse_utils
from data_warehouse_utils.dataloader import DataLoader

from psycopg2 import sql
from data_warehouse_utils import config, connect
from data_warehouse_utils.utils import get_db_config

import sqlalchemy
from sqlalchemy.engine import Engine
import pandas as pd
import os

db_config = get_db_config()
engine = connect.get_db_engine(db_name=db_config['db_name'], ip=db_config['ip'],
                               username = db_config['username'], 
                               password = db_config['password'])
# pd.read_sql_query(
#   sql = sql.SQL('SELECT COUNT(*) FROM processed.single_timestamp'),
#   con=engine.raw_connection()
# )
convert_ST = True

names = [
  'episodes', 'admissions',
  'medications', 'intubations', 'range_measurements',
  'diagnoses', 'parameters', 'comorbidities', 'outcomes'
]

for tbl_name in names:
  df = pd.read_sql_query(
    sql = sql.SQL(f'SELECT * FROM processed.{tbl_name}'),
    con=engine.raw_connection()
  )
  df.to_parquet(tbl_name + '.parquet')

# df = pd.read_sql_query(
#     sql = sql.SQL(f'SELECT * FROM processed.episodes'),
#     con=engine.raw_connection()
# )

if convert_ST:
  
  hospitals = [
    'stantonius', 'erasmus', 'umcu', 'cwz', 'radboud', 'catharina',
       'olvg', 'jeroenbosch', 'amc', 'albertschweitzer', 'franciscus',
       'spaarnegasthuis', 'zuyderland', 'martini', 'haga', 'noordwest',
       'maasstad', 'amphia', 'viecuri', 'vumc', 'etz', 'bovenij', 'rdgg',
       'laurentius', 'geldersevallei', 'ikazia', 'ijsselland', 'diakhuis',
       'zgt', 'wza'
  ]
  
  if not os.path.exists("single_timestamp"):
    os.mkdir("single_timestamp")
  
  curr_size = 1
  idx = 0
  tbl_name = 'single_timestamp'
  
  for hosp in hospitals:
    
    df = pd.read_sql_query(
      sql = sql.SQL(f'SELECT * FROM processed.{tbl_name} WHERE hospital = \'{hosp}\' '),
      con=engine.raw_connection()
    )
    
    idx = idx + 1
    curr_size = df.size
    if df.size > 0:
      df.to_parquet('single_timestamp/' + str(idx) + '.parquet')


# dl = DataLoader()
# pts = dl.get_episodes()

# df = pd.read_sql_query(
#   sql = sql.SQL(f'SELECT * FROM processed.single_timestamp WHERE hospital = \'stantonius\' '),
#   con=engine.raw_connection()
# )
