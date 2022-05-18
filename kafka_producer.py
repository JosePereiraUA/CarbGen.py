from __future__ import print_function
from kafka import KafkaProducer
from kafka.errors import KafkaError
import sys
import json


topic = 'test'
server = '192.168.179.53:9092'


producer = KafkaProducer(
             bootstrap_servers=server,
             value_serializer=lambda v: json.dumps(v).encode('utf-8')
           )

msg =1837.7 if len(sys.argv)==1 else sys.argv[1]

future = producer.send(topic, msg)
#producer.flush()
#result = future.get(timeout=60)
#print(result)
producer.close()

