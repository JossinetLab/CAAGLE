var express = require('express'),
autocomplete=require('autocomplete');


var app = express();

app.configure(function () {
 app.use(express.logger('dev')); /* 'default', 'short', 'tiny', 'dev' */
 app.use(express.bodyParser());
});


app.get('/autocomplete/:search',autocomplete.find);

app.listen(6000);
console.log('Listening on port 3000...');
