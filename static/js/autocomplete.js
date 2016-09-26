var mongo = require('mongodb');

var Server = mongo.Server,
        Db = mongo.Db,
        BSON = mongo.BSONPure;

var server = new Server('localhost', 27017, {auto_reconnect: true});
db = new Db('redsocial', server);

db.open(function(err, db) {
    if (!err) {
        console.log("Connected to 'mydb' database");
        db.collection('publication', {strict: true}, function(err, collection) {
            if (err) {
                console.log("error");
            }
        });
    }
});


exports.find = function(req, res) {
var b=req.params.search;
db.collection('publication', function(err, collection) {
      collection.find({type:'pub',content: new RegExp(b,'i')}).limit(5).toArray(function(err, items) {
                res.jsonp(items);
            });
        });
};
