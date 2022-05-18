from app import db


class Job(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    completed = db.Column(db.Boolean, default=False)
    result = db.Column(db.BLOB, nullable=True)

    def __str__(self):
        return '<Job {}:{}>'.format(self.id, 'done' if self.completed else 'pending')
    
    def __repr__(self):
        return str(self) 
    
class Statistics(db.Model):
    country = db.Column(db.String(120), primary_key=True, nullable=False)
    count = db.Column(db.Integer, default=0)

    def __init__(self, country, count=0):
        db.Model.__init__(self)
        self.country = country
        self.count = count

    def __str__(self):
        return '<Statistics {}:{}>'.format(self.country, self.count)
    
    def __repr__(self):
        return str(self) 