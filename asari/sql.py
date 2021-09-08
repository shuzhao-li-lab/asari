'''

Use SQLite3 as backend database to store all sample-related data.

This relieves the memory requirement and provides persistence to stop/resume operations.

Multiple processes can visit the same database, thus enabling parallel processing.


Design
------

One Experiment (project) per database.

Each sample links to two tables: mass_traces and peaks.






Use pandas for DB operations (ORM).

Use JOIN etc for data operation.






For ref DBs, flat the relationship list and use the 1st only in SQL table.
[['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]],





'''


import sqlite3 as sql
import pandas as pd





