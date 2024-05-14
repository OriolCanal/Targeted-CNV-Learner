from bson import ObjectId


def get_timestamp(Run):
    # If annotation is a DBRef object, return None
    if isinstance(Run, ObjectId):
        return None
    
    # If annotation is an object with a 'timestamp' attribute
    if hasattr(Run, "timestamp"):
        return Run.timestamp
    


    